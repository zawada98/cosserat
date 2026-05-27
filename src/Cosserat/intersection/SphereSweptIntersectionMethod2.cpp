/******************************************************************************
 * Cosserat Plugin for SOFA Framework                                         *
 *                                                                            *
 * SphereSweptIntersectionMethod.cpp                                          *
 *                                                                            *
 * See SphereSweptIntersectionMethod.h for full documentation.               *
 ******************************************************************************/
#include <fstream>
#include <cstdlib>

#include "SphereSweptIntersectionMethod.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/accessor.h>
#include <algorithm>

namespace Cosserat {

    namespace
    {
        bool traceZeroContactNormals()
        {
            static const bool enabled = std::getenv("COSSERAT_TRACE_ZERO_CONTACT_NORMALS") != nullptr;
            return enabled;
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Constructor
    // ─────────────────────────────────────────────────────────────────────────────
    SphereSweptIntersectionMethod::SphereSweptIntersectionMethod()
        : Inherit1()
        //
        // ── Inputs ────────────────────────────────────────────────────────────
        , d_beam1Frames(initData(&d_beam1Frames, VecRigidCoord{},
            "beam1Frames",
            "Rigid3d frames of Beam 1 (FramesMO output of DiscreteCosseratMapping)."))
        , d_beam2Frames(initData(&d_beam2Frames, VecRigidCoord{},
            "beam2Frames",
            "Rigid3d frames of Beam 2 (FramesMO output of DiscreteCosseratMapping)."))
        , d_beam1Velocities(initData(&d_beam1Velocities, VecRigidDeriv{},
            "beam1Velocities",
            "Rigid3d frame velocities of Beam 1 (FramesMO.velocity). "
            "Link via beam1Velocities=beam1_MO.getLinkPath()+'.velocity'."))
        , d_beam2Velocities(initData(&d_beam2Velocities, VecRigidDeriv{},
            "beam2Velocities",
            "Rigid3d frame velocities of Beam 2 (FramesMO.velocity). "
            "Same rules as beam1Velocities."))
        , d_radius1(initData(&d_radius1, Real(0.1),
            "radius1",
            "Outer cross-section radius of Beam 1."))
        , d_radius2(initData(&d_radius2, Real(0.1),
            "radius2",
            "Outer cross-section radius of Beam 2."))
        , d_innerRadius1(initData(&d_innerRadius1, Real(0.0),
            "innerRadius1",
            "Inner radius of Beam 1. 0 = solid beam (external contact). "
            ">0 = hollow tube (CTR nesting). "
            "Beam 1 is the CTR outer tube when radius1 > radius2."))
        , d_innerRadius2(initData(&d_innerRadius2, Real(0.0),
            "innerRadius2",
            "Inner radius of Beam 2. 0 = solid beam (external contact). "
            ">0 = hollow tube (CTR nesting). "
            "Beam 2 is the CTR outer tube when radius2 > radius1."))
        , d_contactConfiguration(initData(&d_contactConfiguration,
            sofa::helper::OptionsGroup({ "external", "nested" }),
            "contactConfiguration",
            "Contact geometry mode. "
            "'external': always side-by-side; gap = dist-(r1+r2). "
            "'nested': always coaxial (CTR); gap = ri_outer-r_outer_inner-dist. "
            "  One-time init check: ri_outer > r_outer_inner must hold. "
            "  defaultNormal MUST be set to a vector perpendicular to the beam axis."))
        , d_broadPhaseMarginFactor(initData(&d_broadPhaseMarginFactor, Real(1.5),
            "broadPhaseMarginFactor",
            "Multiplier on each beam's contact-relevant radius for bounding-sphere broad-phase. "
            "External: R_i = halfLen_i + factor*r_outer_i. "
            "Nested (outer tube): R_i = halfLen_i + factor*ri_outer_i (bore). "
            "Nested (inner tube): R_i = halfLen_i + factor*r_outer_i. "
            "Default 1.5. For nested CTR beams all segments pass (correct)."))
        , d_defaultNormal(initData(&d_defaultNormal, Vec3(0,0,0),
            "defaultNormal",
            "Fallback contact normal for degenerate cases (coincident + parallel beams "
            "with no prior valid normal). Must be set by the user. "
            "If left as (0,0,0) the contact pair is skipped with an error. "
            "REQUIRED for nested CTR scenes (centrelines nearly coincident by design). "
            "Use a vector perpendicular to the nominal beam axis (e.g. '0 1 0')."))
        , d_cachedNormalMaxAxialProjection(initData(&d_cachedNormalMaxAxialProjection, 
            Real(0.17),
            "cachedNormalMaxAxialProjection",
            "Tolerance on |n_cached . tangent| for cached-normal perpendicularity test. "
            "A cached normal is rejected if its projection on either segment tangent "
            "exceeds this value. Default 0.17 ~ sin(10 deg)."))
        // ── Outputs ───────────────────────────────────────────────────────────
        , d_curvilinearParams(initData(&d_curvilinearParams,
            "curvilinearParams",
            "Normalised curvilinear parameters {s1*, s2*} per contact pair. "
            "Always refers to original beam numbering. "))
    {
        sofa::helper::OptionsGroup cfgOptions({ "external", "nested" }); 
        cfgOptions.setSelectedItem(0u); // default: "external"
        d_contactConfiguration.setValue(cfgOptions);
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  SOFA lifecycle
    // ─────────────────────────────────────────────────────────────────────────────
    void SphereSweptIntersectionMethod::init()
    {
        addInput(&d_beam1Frames);
        addInput(&d_beam2Frames);
        addInput(&d_radius1);
        addInput(&d_radius2);
        addInput(&d_innerRadius1);
        addInput(&d_innerRadius2);
        addInput(&d_contactConfiguration);
        addInput(&d_broadPhaseMarginFactor);
        addInput(&d_defaultNormal);
        addInput(&d_beam1Velocities);  
        addInput(&d_beam2Velocities); 

        addOutput(&d_curvilinearParams);
        
        Inherit1::init();
       	validateParameters();
        setDirtyValue();
    }

    void SphereSweptIntersectionMethod::reinit()
    {
        if (!validateParameters())
            return;
        update();
    }
    
// ─────────────────────────────────────────────────────────────────────────────
//  validateParameters
//
//  Runs every parameter-consistency check. Sets d_componentState to Invalid on
//  failure and returns false. Clears m_lastValidNormal because cached normals
//  keyed by (i,j) section-pair indices may be stale after topology changes.
//  Called by init() (after addInput/addOutput) and reinit().
// ─────────────────────────────────────────────────────────────────────────────
bool SphereSweptIntersectionMethod::validateParameters()
{
    // Reset state; failure branches set it to Invalid.
    d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);

    // Drop cache — geometry/topology may have changed since last validation.
    m_lastValidNormal.clear();

    // ── defaultNormal check ──────────────────────────────────────────────
    if (d_defaultNormal.getValue().norm() < s_eps)
    {
        const std::string cfg =
            d_contactConfiguration.getValue().getSelectedItem();
        if (cfg == "nested")
        {
            msg_error() << "contactConfiguration='nested' requires 'defaultNormal' "
                           "to be a non-zero vector perpendicular to the beam axis.";
            d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
            return false;
        }
        else
        {
            msg_warning() << "'defaultNormal' is zero. Used only as last-resort "
                             "fallback for coincident+parallel beams.";
        }
    }

    if (!d_beam1Frames.getParent() || !d_beam2Frames.getParent())
    {
        msg_error() << "beam1Frames and beam2Frames must be linked to a MechanicalObject.";
        d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }
    if (d_radius1.getValue() <= Real(0) || d_radius2.getValue() <= Real(0))
    {
        msg_error() << "radius1 and radius2 must be strictly positive.";
        d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }
    if (!d_beam1Velocities.getParent())
    {
        msg_error() << "'beam1Velocities' is MANDATORY. Link via "
                       "beam1Velocities=<beam1_framesMO>.getLinkPath()+'.velocity'.";
        d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }
    if (!d_beam2Velocities.getParent())
    {
        msg_error() << "'beam2Velocities' is MANDATORY. Link via "
                       "beam2Velocities=<beam2_framesMO>.getLinkPath()+'.velocity'.";
        d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }

    sofa::helper::ReadAccessor<sofa::Data<VecRigidCoord>> f1 = d_beam1Frames;
    sofa::helper::ReadAccessor<sofa::Data<VecRigidCoord>> f2 = d_beam2Frames;
    sofa::helper::ReadAccessor<sofa::Data<VecRigidDeriv>> v1 = d_beam1Velocities;
    sofa::helper::ReadAccessor<sofa::Data<VecRigidDeriv>> v2 = d_beam2Velocities;

    if (v1.size() != f1.size() || v2.size() != f2.size())
    {
        msg_error() << "Velocity array size does not match frame array size "
                       "(beam1: " << v1.size() << " vs " << f1.size()
                    << ", beam2: " << v2.size() << " vs " << f2.size() << ").";
        d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }

    if (d_innerRadius1.getValue() < Real(0) || d_innerRadius2.getValue() < Real(0))
    {
        msg_error() << "innerRadius1 and innerRadius2 must be non-negative.";
        d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }
    if (d_innerRadius1.getValue() > Real(0) &&
        d_innerRadius1.getValue() >= d_radius1.getValue())
    {
        msg_error() << "innerRadius1 (" << d_innerRadius1.getValue()
                    << ") must be strictly less than radius1 ("
                    << d_radius1.getValue() << ").";
        d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }
    if (d_innerRadius2.getValue() > Real(0) &&
        d_innerRadius2.getValue() >= d_radius2.getValue())
    {
        msg_error() << "innerRadius2 (" << d_innerRadius2.getValue()
                    << ") must be strictly less than radius2 ("
                    << d_radius2.getValue() << ").";
        d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }
        // ── Broad-phase margin ──────────────────────────────────────────────
    const Real bpf = d_broadPhaseMarginFactor.getValue();
    if (bpf <= Real(0))
    {
        msg_error() << "broadPhaseMarginFactor (" << bpf
                    << ") must be strictly positive. "
                       "Typical value: 1.5.";
        d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }
    if (bpf < Real(1))
    {
        msg_warning() << "broadPhaseMarginFactor (" << bpf
                      << ") < 1 shrinks the bounding spheres below the "
                         "true tube radii and will miss real contacts. "
                         "Recommended: >= 1 (default 1.5).";
    }

    // ── Cached-normal perpendicularity tolerance ────────────────────────
    // d_cachedNormalMaxAxialProjection is compared against |n_cached . tangent|
    // where both operands are unit vectors, so the meaningful range is [0, 1].
    const Real axTol = d_cachedNormalMaxAxialProjection.getValue();
    if (axTol < Real(0) || axTol > Real(1))
    {
        msg_error() << "cachedNormalMaxAxialProjection (" << axTol
                    << ") must lie in [0, 1] (it is compared to "
                       "|n . tangent| of unit vectors). "
                       "Default: 0.17 = sin(10 deg).";
        d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }
    if (axTol > Real(0.5))   // sin(30 deg)
    {
        msg_warning() << "cachedNormalMaxAxialProjection (" << axTol
                      << ") > 0.5 (= sin 30 deg) accepts cached normals "
                         "that are nearly axial. Smoothness gate is weak.";
    }


    // ── Nested geometry consistency ──────────────────────────────────────
    const std::string cfg = d_contactConfiguration.getValue().getSelectedItem();
    if (cfg == "nested")
    {
        const Real r1  = d_radius1.getValue(),      r2  = d_radius2.getValue();
        const Real ri1 = d_innerRadius1.getValue(), ri2 = d_innerRadius2.getValue();
        const bool beam1IsOuter  = (r1 >= r2);
        const Real r_bore        = beam1IsOuter ? ri1 : ri2;
        const Real r_outer_inner = beam1IsOuter ? r2  : r1;
        const int  outerBeamId   = beam1IsOuter ? 1   : 2;

        if (r_bore <= Real(0))
        {
            msg_error() << "contactConfiguration='nested': Beam " << outerBeamId
                        << " is the outer tube but its innerRadius is 0. "
                           "Set innerRadius" << outerBeamId << " > 0.";
            d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
            return false;
        }
        if (r_outer_inner >= r_bore)
        {
            msg_error() << "contactConfiguration='nested': inner tube outer radius ("
                        << r_outer_inner << ") >= outer tube bore (" << r_bore
                        << "). Inner tube cannot fit at construction.";
            d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
            return false;
        }
    }

    return true;
}

    // ─────────────────────────────────────────────────────────────────────────────
    //  doUpdate  – main entry point
    // ─────────────────────────────────────────────────────────────────────────────
    std::ofstream& ssimLog() {
        // Opens once. Path is overridable via the COSSERAT_LOG_DIR env var,
        // otherwise drops next to wherever runSofa was launched from.
        static std::ofstream f = [] {
            const char* dir = std::getenv("COSSERAT_LOG_DIR");
            std::string path = (dir ? std::string(dir) + "/" : std::string("")) + "ssim_log.txt";
            return std::ofstream(path, std::ios::out | std::ios::trunc);
        }();
        return f;
    }
    
    void SphereSweptIntersectionMethod::doUpdate()
    {   
        // ── Write-accessors for outputs ──────────────────────────────────────────
        sofa::helper::WriteOnlyAccessor<sofa::Data<sofa::type::vector<sofa::type::Vec2d>>>
            outParams = d_curvilinearParams;
        
        m_distances.clear();         
        m_centerlinePoints1.clear();  
        m_centerlinePoints2.clear(); 
        m_surfacePoints1.clear();     
        m_surfacePoints2.clear();     
        m_contactSectionIds.clear();    
        m_contactNormals.clear();
        m_contactTangents1.clear();
        m_contactTangents2.clear();

        outParams.clear();
        
        if (this->d_componentState.getValue() ==
        sofa::core::objectmodel::ComponentState::Invalid)
        {
            return;
        }
        
        // ── Read inputs ──────────────────────────────────────────────────────────
        sofa::helper::ReadAccessor<sofa::Data<VecRigidCoord>> frames1 = d_beam1Frames;
        sofa::helper::ReadAccessor<sofa::Data<VecRigidCoord>> frames2 = d_beam2Frames;
        
        sofa::helper::ReadAccessor<sofa::Data<VecRigidDeriv>> vels1 = d_beam1Velocities;
        sofa::helper::ReadAccessor<sofa::Data<VecRigidDeriv>> vels2 = d_beam2Velocities;

        // dt is needed to convert tangential velocity → tangential displacement.
        const Real dt_sim = static_cast<Real>(getContext()->getDt());

        const Real r1     = d_radius1.getValue();
        const Real r2     = d_radius2.getValue();
        const Real ri1    = d_innerRadius1.getValue();
        const Real ri2    = d_innerRadius2.getValue();
        const Real margin = d_broadPhaseMarginFactor.getValue();

        const int N1 = static_cast<int>(frames1.size());
        const int N2 = static_cast<int>(frames2.size());

        if (N1 < 2 || N2 < 2)
        {
            msg_warning() << "Need at least 2 frames per beam. "
                "Beam1 has " << N1 << ", Beam2 has " << N2 << ".";
            return;
        }

        // ── Contact configuration ────────────────────────────────────────────────
        //
        // "external": always side-by-side, regardless of radii.
        //    gap = dist − (r1 + r2).
        //
        // "nested": always coaxial (CTR). The beam with the larger outer radius is
        //    the outer tube (determined once here, never changes per segment).
        //    gap = ri_outer − r_outer_inner − dist.
        //    Validity already verified in init().
        //
        const std::string contactCfg =
            d_contactConfiguration.getValue().getSelectedItem();

        const bool isNested = (contactCfg == "nested");

        const bool beam1IsOuter  = (r1 >= r2);
        const Real r_bore        = beam1IsOuter ? ri1 : ri2; // inner radius of outer tube (bore)
        const Real r_outer_inner = beam1IsOuter ? r2  : r1;  // outer radius of inner tube

        // Gap lambda: takes the centreline distance for the current segment pair.
        auto computeGap = [&](Real dist) -> Real
        {
            if (isNested)
                return r_bore - r_outer_inner - dist;
            return dist - (r1 + r2); // "external"
        };
       

        // The finer beam (more segments / nodes) drives the outer loop to avoid
        // missing contacts. When swapped (Beam 2 is finer), the internal segment
        // roles are swapped but all outputs are re-labelled to the original beam
        // numbering before being written.
        const bool useSwapped = (N2 > N1);

        const VecRigidCoord& finerFrames   = useSwapped ? frames2.ref() : frames1.ref();
        const VecRigidCoord& coarserFrames = useSwapped ? frames1.ref() : frames2.ref();
        const Real r_finer   = useSwapped ? r2 : r1; // outer radius of finer-mesh beam
        const Real r_coarser = useSwapped ? r1 : r2; // outer radius of coarser-mesh beam
        const Real ri_finer  = useSwapped ? ri2 : ri1; // inner radius of finer-mesh beam
        const Real ri_coarser= useSwapped ? ri1 : ri2; // inner radius of coarser-mesh beam

        const int N_finer   = static_cast<int>(finerFrames.size());
        const int N_coarser = static_cast<int>(coarserFrames.size());

        // determine whether the finer beam is the outer tube.
        // This is independent of the finer/coarser distinction (mesh resolution
        // does not necessarily correlate with tube nesting).
        const bool finerIsOuter = useSwapped ? !beam1IsOuter : beam1IsOuter;

        // contact-relevant broad-phase radii.
        // For the outer tube (nested mode), the bore (ri) is the contact surface;
        // using the outer radius would over-inflate the bounding sphere.
        // For external mode (or the inner tube in nested mode), use the outer radius.
        const Real r_bp_query     = (isNested && finerIsOuter)  ? ri_finer   : r_finer;
        const Real r_bp_candidate = (isNested && !finerIsOuter) ? ri_coarser : r_coarser;

        // compute mapToSurface radii based on contact configuration
        // and which beam is outer/inner.
        // mapToSurface formula (using raw Beam-1→Beam-2 normal nHat):
        //   psurf1 = pint1 + r_map1 * nHat
        //   psurf2 = pint2 +/- r_map2 * nHat (+ if nested, - otherwise)
        //
        // External: outer surfaces face each other → r_map1=r1, r_map2=r2.
        // Nested (beam1=outer): psurf1 on bore (ri1), psurf2 on outer wall (r2).
        // Nested (beam1=inner): psurf1 on outer wall (r1), psurf2 on bore (ri2).
        const Real r_map1 = !isNested       ?  r1      // external: +r1
                   :  beam1IsOuter   ?  ri1     // nested, B1=outer: +ri1 (bore, on +n̂ side)
                   :                    -r1;    // nested, B1=inner: −r1  (wall, on −n̂ side)

        const Real r_map2 = !isNested       ?  r2      // external: +r2  (produces −r2·n̂)
                          :  beam1IsOuter   ?  -r2     // nested, B2=inner: −r2  (so psurf2 = Pb + r2·n̂)
                          :                    ri2;    // nested, B2=outer: +ri2 (so psurf2 = Pb − ri2·n̂)
        
        for (int i = 0; i < N_finer - 1; ++i)
        {
            const Vec3 p0 = finerFrames[i].getCenter();
            const Vec3 p1 = finerFrames[i + 1].getCenter();

            const Vec3 seg_tangent_finer = [&]() -> Vec3 {
                const Vec3 d = p1 - p0;
                const Real dn = d.norm();
                return (dn > s_eps) ? d / dn
                    : finerFrames[i].getOrientation().rotate(Vec3(Real(1), Real(0), Real(0)));
            }();

            const Vec3 mid_i     = (p0 + p1) * Real(0.5);
            const Real halfLen_i = (p1 - p0).norm() * Real(0.5);

            const auto candidates = candidateSegments(
                mid_i, halfLen_i, r_bp_query, coarserFrames, r_bp_candidate, margin);

            // Best-candidate selection differs by contact configuration:                  
            //   external: smallest centreline distance wins (Ericson §5.1.9).             
            //   nested  : axial-overlap filter + LARGEST radial separation over the      
            //             overlap wins. Candidates with no axial overlap on this T1      
            //             segment are rejected — they are in axial free space relative   
            //             to this T1 segment and would produce spurious gap values.       
            Real bestDist = isNested                                                      
                            ? std::numeric_limits<Real>::lowest()                         
                            : std::numeric_limits<Real>::max();                     
            Real best_s_o = Real(0), best_s_i = Real(0);
            Vec3 best_cp_o, best_cp_i;
            int  best_j   = -1;
            bool foundAny = false;

            for (const int j : candidates)
            {
                const Vec3 q0 = coarserFrames[j].getCenter();
                const Vec3 q1 = coarserFrames[j + 1].getCenter();

                Real s_o = Real(0), s_i = Real(0);
                Vec3 cp_o, cp_i;

                if (!isNested)                                                             
                {                                                                         
                    // ── External: standard segment-to-segment min distance ────────────
                    segmentToSegment(p0, p1, q0, q1, s_o, s_i, cp_o, cp_i);

                    Real centrelineDist = (cp_o - cp_i).norm();
                    if (centrelineDist < s_eps)
                    {
                        msg_error() << "External mode, centerlines coincide at segment pair ("  
                                    << (useSwapped ? j : i) << ","                                       // (msg_error path is now only reachable in external mode;
                                    << (useSwapped ? i : j)                                              //  the original `if (!isNested)` guard around the message
                                    << "). This is a real interpenetration.";                           //  is therefore redundant and removed.)
                        cp_o = (p0 + p1) * Real(0.5);
                        nodeToSegment(cp_o, q0, q1, s_i, cp_i);
                        s_o  = Real(0.5);
                        centrelineDist = (cp_o - cp_i).norm();
                    }

                    if (centrelineDist < bestDist)
                    {
                        bestDist  = centrelineDist;
                        best_s_o  = s_o; best_s_i = s_i;
                        best_cp_o = cp_o; best_cp_i = cp_i;
                        best_j    = j;
                        foundAny  = true;
                    }
                }                                                                          
                else                                                                      
                {                                                                          
                    // ── Nested: axial-overlap filter + max-radial-in-overlap ──────────  
                    Real radialDist;                                                       
                    if (!axialOverlapMaxRadial(p0, p1, q0, q1,                            
                                               s_o, s_i, cp_o, cp_i, radialDist))          
                    {                                                                      
                        continue;  // T2 seg j axially outside this T1 segment.            
                    }                                                                      
                                                                                  
                    if (radialDist > bestDist)                                           
                    {                                                                   
                        bestDist  = radialDist;                                            
                        best_s_o  = s_o; best_s_i = s_i;                                
                        best_cp_o = cp_o; best_cp_i = cp_i;                            
                        best_j    = j;                                                   
                        foundAny  = true;                                                
                    }                                                                     
                }                                                                         
            }
            
            if (!foundAny) continue;
            const Real bestGap = computeGap(bestDist);

            // Re-label to original beam-1/beam-2 numbering.
            const int  idx1   = useSwapped ? best_j : i;
            const int  idx2   = useSwapped ? i      : best_j;
            const Vec3 pint1  = useSwapped ? best_cp_i : best_cp_o;
            const Vec3 pint2  = useSwapped ? best_cp_o : best_cp_i;
            const Real s1_out = useSwapped ? best_s_i  : best_s_o;
            const Real s2_out = useSwapped ? best_s_o  : best_s_i;

            const Vec3 q0_best = coarserFrames[best_j].getCenter();
            const Vec3 q1_best = coarserFrames[best_j + 1].getCenter();
            const Vec3 seg_tangent_coarser = [&]() -> Vec3 {
                const Vec3 d = q1_best - q0_best;
                const Real dn = d.norm();
                return (dn > s_eps) ? d / dn
                    : coarserFrames[best_j].getOrientation().rotate(Vec3(Real(1), Real(0), Real(0)));
            }();

            // Reassign to original beam-1/beam-2 tangents after the swap.
            const Vec3& t1 = useSwapped ? seg_tangent_coarser : seg_tangent_finer;
            const Vec3& t2 = useSwapped ? seg_tangent_finer   : seg_tangent_coarser;

            // Compute raw contact normal (always Beam-1 → Beam-2 direction).
            const Vec3 nHat_raw = computeContactNormal(
                pint1, pint2, idx1, idx2, t1, t2,
                frames1.ref()[idx1], frames1.ref()[idx1 + 1], s1_out);

            // for nested mode, enforce the outer→inner sign convention
            // for the output normal, regardless of which beam is Beam-1 or Beam-2.
            // nHat_raw = Beam-1 → Beam-2.
            // If Beam-1 is the inner tube (!beam1IsOuter), nHat_raw points inner→outer;
            // flip it so the reported normal always points outer→inner.
            // mapToSurface uses nHat_raw (Beam-1→Beam-2) internally — the r_map1/r_map2
            // values are already chosen to be consistent with this convention.
            // Previously: no sign flip; output normal was Beam-1→Beam-2 regardless of mode.
            const Vec3 nHat_out = (isNested && !beam1IsOuter) ? -nHat_raw : nHat_raw;
            m_contactNormals.push_back(nHat_out);

            // mapToSurface uses correct contact-relevant radii (r_map1, r_map2)
            // instead of always using the outer radii (r1, r2).
            // For nested (beam1=outer): r_map1=ri1 (bore), r_map2=r2 (outer wall of inner).
            // For nested (beam1=inner): r_map1=r1 (outer wall of inner), r_map2=ri2 (bore).
            // For external: r_map1=r1, r_map2=r2 (unchanged from before).
            // Previously: mapToSurface(pint1, pint2, nHat, r1, r2, psurf1, psurf2);
            Vec3 psurf1, psurf2;
            mapToSurface(pint1, pint2, nHat_raw, r_map1, r_map2, psurf1, psurf2);

            // ── Contact-plane tangent basis ───────────────────────────────────
            Vec3 t1_contact, t2_contact;
            computeContactFrame(t1 /*tau1*/, t2 /*tau2*/, nHat_out, t1_contact, t2_contact);

            m_contactTangents1.push_back(t1_contact);
            m_contactTangents2.push_back(t2_contact);

            Real delta_t1 = Real(0);                                             
            Real delta_t2 = Real(0);                                        

            
            const int next1 = std::min(idx1 + 1, (int)vels1.size() - 1);
            const int next2 = std::min(idx2 + 1, (int)vels2.size() - 1);

            // Frame centres (needed for moment arms r = Pa/Pb − p_frame)  
            const Vec3 pf1_i    = frames1.ref()[idx1].getCenter();
            const Vec3 pf1_next = frames1.ref()[next1].getCenter();
            const Vec3 pf2_i    = frames2.ref()[idx2].getCenter();
            const Vec3 pf2_next = frames2.ref()[next2].getCenter();

            // Angular velocities of the bounding frames             
            const Vec3 omega1_i    = vels1[idx1].getVOrientation();
            const Vec3 omega1_next = vels1[next1].getVOrientation();
            const Vec3 omega2_i    = vels2[idx2].getVOrientation();
            const Vec3 omega2_next = vels2[next2].getVOrientation();

            // Moment arms: r = contact_point − frame_centre               
            const Vec3 r1_i    = psurf1 - pf1_i;
            const Vec3 r1_next = psurf1 - pf1_next;
            const Vec3 r2_i    = psurf2 - pf2_i;
            const Vec3 r2_next = psurf2 - pf2_next;

            // ω × r  
            auto cross = [](const Vec3& a, const Vec3& b) -> Vec3 {
                return Vec3(a[1]*b[2] - a[2]*b[1],
                            a[2]*b[0] - a[0]*b[2],
                            a[0]*b[1] - a[1]*b[0]);
            };

            // Full rigid-body velocity at the contact point:               
            // v_contact = v_centre + ω × r, interpolated across the two    
            // bounding frames at the contact curvilinear parameter.        
            // Previously: v_Pa = slerp(v_i, v_{i+1}, s) — angular term omitted 
            const Vec3 v_Pa =
                (vels1[idx1].getVCenter()  + cross(omega1_i,    r1_i))    * (Real(1) - s1_out)
              + (vels1[next1].getVCenter() + cross(omega1_next, r1_next)) * s1_out;  

            const Vec3 v_Pb =
                (vels2[idx2].getVCenter()  + cross(omega2_i,    r2_i))    * (Real(1) - s2_out)
              + (vels2[next2].getVCenter() + cross(omega2_next, r2_next)) * s2_out; 

            const Vec3 v_rel = v_Pb - v_Pa;
            delta_t1 = (v_rel * t1_contact) * dt_sim;                       
            delta_t2 = (v_rel * t2_contact) * dt_sim;
            
            m_distances.push_back(Vec3(bestGap, delta_t1, delta_t2));
            m_centerlinePoints1.push_back(pint1);                       
            m_centerlinePoints2.push_back(pint2);                      
            m_surfacePoints1.push_back(psurf1);                        
            m_surfacePoints2.push_back(psurf2);                       
            m_contactSectionIds.push_back({ idx1, idx2 });       
            
            outParams.push_back({ s1_out, s2_out });
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  computeContactNormal
    //
    //  Returns the RAW Beam-1 → Beam-2 contact normal.
    //  The nested outer→inner sign flip is applied by doUpdate() on the output,
    //  NOT here, so that m_lastValidNormal caches the raw direction and the cache
    //  remains consistent regardless of which beam is outer/inner in the scene.
    //
    //  Zero-normal fallback order:
    //    1. Last valid normal for pair (i,j) from m_lastValidNormal.
    //    2. Cross product of segment tangents t̂1 × t̂2.
    //    3. User-supplied defaultNormal rotated by SLERP frame at contact point.
    // ─────────────────────────────────────────────────────────────────────────────
    SphereSweptIntersectionMethod::Vec3
SphereSweptIntersectionMethod::computeContactNormal(
    const Vec3&       pint1,
    const Vec3&       pint2,
    int               i,
    int               j,
    const Vec3&       seg1Tangent,
    const Vec3&       seg2Tangent,
    const RigidCoord& frameA,
    const RigidCoord& frameB,
    Real              s1)
{
    // Current SLERP'd Beam-1 orientation at the contact point — used both          
    // for caching on the write path and for the smoothness test on the read path.
    // Previously: computed only in fallback 2c.
    const auto qSlerpNow = frameA.getOrientation().slerp(frameB.getOrientation(), s1);

    // ── 1. Standard case: non-degenerate centreline separation ───────────────
    Vec3       n     = pint2 - pint1;
    const Real dnorm = n.norm();

    if (dnorm >= s_eps)
    {
        n /= dnorm;
        m_lastValidNormal[{i, j}] = {n, qSlerpNow};                                
        return n;
    }

    // ── 2. Zero-normal fallback ──────────────────────────────────────────────

    // 2a. Reuse last valid normal for this (i,j) pair, ONLY if it passes         
    //     two validity tests:
    //       (a) perpendicularity: |n_cached · t̂₁| and |n_cached · t̂₂| small
    //       (b) smoothness     : Beam-1 SLERP'd frame rotated by at most
    //                            no snapthrough (180 rotation) since cache time.
    //     If either test fails, we DO NOT use the cached value but we DO NOT
    //     erase it — it may become usable again in a later timestep.
    //     Previously: cached normal was returned unconditionally.
    auto it = m_lastValidNormal.find({i, j});
    if (it != m_lastValidNormal.end())
    {
        const Vec3& n_cached = it->second.normal;
        const auto& q_cached = it->second.qSlerp;

        // Test (a): perpendicularity to both segment tangents.
        const Real tolPerp = d_cachedNormalMaxAxialProjection.getValue();
        const Real proj1   = std::abs(n_cached * seg1Tangent); // Vec3::operator* is dot
        const Real proj2   = std::abs(n_cached * seg2Tangent);
        const bool passPerp = (proj1 <= tolPerp) && (proj2 <= tolPerp);

        // Test (b): smoothness of SLERP'd contact-point frame.
        // Relative quat q_rel = q_now * q_cached^{-1}; angle = 2*acos(|w|).
        const auto qRel    = qSlerpNow * q_cached.inverse();
        Real       wAbs    = std::abs(qRel[3]); // scalar component; convention [0..2]=xyz, [3]=w
        if (wAbs > Real(1)) wAbs = Real(1);     // clamp for acos numerical safety
        const Real rotAngle  = Real(2) * std::acos(wAbs);
        const Real tolRot    = Real(3.14159265358979323846/ 18.0); // just verifies that no snap through occured
        const bool passSmooth = (rotAngle <= tolRot);

        if (passPerp && passSmooth)
        {
            if (traceZeroContactNormals())
            {
                msg_warning() << "Contact pair (" << i << "," << j
                              << "): zero contact normal (coincident centrelines). "
                                 "Reusing last valid normal (validated: |n.t1|=" << proj1
                              << ", |n.t2|=" << proj2
                              << ", frame rot=" << rotAngle << " rad).";
            }
            return n_cached;
        }

        // Cache not trusted for this step — fall through without erasing.
        if (traceZeroContactNormals())
        {
            msg_warning() << "Contact pair (" << i << "," << j
                          << "): zero contact normal, but cached normal rejected "
                              "(|n.t1|=" << proj1 << ", |n.t2|=" << proj2
                          << ", frame rot=" << rotAngle << " rad; tolPerp=" << tolPerp
                          << ", tolRot=" << tolRot
                          << "). Cache preserved; falling back to tangent cross-product.";
        }
    }

    // 2b. Cross product of segment tangents (first occurrence, or cache rejected).
    Vec3 n_cross = Vec3(
        seg1Tangent[1] * seg2Tangent[2] - seg1Tangent[2] * seg2Tangent[1],
        seg1Tangent[2] * seg2Tangent[0] - seg1Tangent[0] * seg2Tangent[2],
        seg1Tangent[0] * seg2Tangent[1] - seg1Tangent[1] * seg2Tangent[0]);
    const Real nc_norm = n_cross.norm();

    if (nc_norm >= s_eps)
    {
        n_cross /= nc_norm;
        if (traceZeroContactNormals())
        {
            msg_warning() << "Contact pair (" << i << "," << j
                          << "): zero contact normal (coincident centrelines). "
                             "Using tangent cross-product fallback.";
        }
        // Only overwrite cache if it was empty; otherwise the rejected entry    
        // might still be useful later once geometry relaxes back.
        // Previously: always overwrote cache.
        if (it == m_lastValidNormal.end())
            m_lastValidNormal[{i, j}] = {n_cross, qSlerpNow};
        return n_cross;
    }

    // 2c. Superimposed parallel beams: rotate user-supplied defaultNormal
    //     by the SLERP frame at the contact point.
    const Vec3 defN = d_defaultNormal.getValue();
    if (defN.norm() < s_eps)
    {
        msg_error() << "Contact pair (" << i << "," << j
                    << "): superimposed parallel beams and defaultNormal is zero. "
                       "Using global Y-axis.";
        return Vec3(Real(0), Real(1), Real(0));
    }
    const Vec3 n_local = qSlerpNow.rotate(defN.normalized());                  
    if (traceZeroContactNormals())
    {
        msg_warning() << "Contact pair (" << i << "," << j
                      << "): superimposed parallel beams, cached normal unavailable/rejected. "
                         "Using user-supplied defaultNormal rotated by SLERP frame.";
    }

    if (it == m_lastValidNormal.end())                                             
        m_lastValidNormal[{i, j}] = {n_local, qSlerpNow};
    return n_local;
}

    // ─────────────────────────────────────────────────────────────────────────────
    //  getContactNormal – accessor for the k-th output slot normal.
    // ─────────────────────────────────────────────────────────────────────────────
    SphereSweptIntersectionMethod::Vec3 SphereSweptIntersectionMethod::getContactNormal(std::size_t k) const
    {   
        ensureUpdated(); 
        if (k >= m_contactNormals.size())
        {
            msg_error() << "getContactNormal: index " << k
                        << " is out of range (current contact count = "
                        << m_contactNormals.size() << "). Aborting.";
            return Vec3(Real(0), Real(0), Real(1));
        }
        return m_contactNormals[k];
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  mapToSurface
    //  Maps centreline contact points to physical surface contact points using the
    //  raw Beam-1→Beam-2 normal nHat and caller-supplied contact-relevant radii:
    //
    //    psurf1 = pint1 + r_surf1 * nHat
    //    psurf2 = pint2 − r_surf2 * nHat
    //
    //  The caller (doUpdate) selects r_surf1/r_surf2 based on the contact mode:
    //    External:            r_surf1=r1,   r_surf2=r2
    //    Nested (beam1=outer):r_surf1=ri1,  r_surf2=r2
    //    Nested (beam1=inner):r_surf1=r1,   r_surf2=ri2
    // ─────────────────────────────────────────────────────────────────────────────
    void SphereSweptIntersectionMethod::mapToSurface(
        const Vec3& pint1,
        const Vec3& pint2,
        const Vec3& nHat,
        Real  r_surf1, Real r_surf2,
        Vec3& psurf1,
        Vec3& psurf2)
    {
        psurf1 = pint1 + nHat * r_surf1;
        psurf2 = pint2 - nHat * r_surf2;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Segment-to-Segment (Ericson 2005, §5.1.9)
    // ─────────────────────────────────────────────────────────────────────────────
    bool SphereSweptIntersectionMethod::segmentToSegment(
        const Vec3& p0, const Vec3& p1,
        const Vec3& q0, const Vec3& q1,
        Real& s1, Real& s2,
        Vec3& cp1, Vec3& cp2)
    {
        const Vec3 d1 = p1 - p0;
        const Vec3 d2 = q1 - q0;
        const Vec3 r  = p0 - q0;

        const Real a = d1.norm2();
        const Real e = d2.norm2();
        const Real f = d2 * r;

        if (a <= s_eps && e <= s_eps)
        {
            s1 = s2 = Real(0);
            cp1 = p0; cp2 = q0;
            return true;
        }

        if (a <= s_eps)
        {
            s1 = Real(0);
            s2 = std::clamp(f / e, Real(0), Real(1));
        }
        else
        {
            const Real c = d1 * r;
            if (e <= s_eps)
            {
                s2 = Real(0);
                s1 = std::clamp(-c / a, Real(0), Real(1));
            }
            else
            {
                const Real b     = d1 * d2;
                const Real denom = a * e - b * b;

                if (denom > s_eps)
                    s1 = std::clamp((b * f - c * e) / denom, Real(0), Real(1));
                else
                    s1 = Real(0);  // parallel segments

                s2 = (b * s1 + f) / e;

                if (s2 < Real(0))
                {
                    s2 = Real(0);
                    s1 = std::clamp(-c / a, Real(0), Real(1));
                }
                else if (s2 > Real(1))
                {
                    s2 = Real(1);
                    s1 = std::clamp((b - c) / a, Real(0), Real(1));
                }
            }
        }

        cp1 = p0 + d1 * s1;
        cp2 = q0 + d2 * s2;
        return true;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  For a LINEAR segment Q(s) = q0 + s*(q1−q0):
    //    s* = clamp( (P−q0)·(q1−q0) / ||q1−q0||² , 0, 1 )
    // ─────────────────────────────────────────────────────────────────────────────
    bool SphereSweptIntersectionMethod::nodeToSegment(
        const Vec3& node,
        const Vec3& q0, const Vec3& q1,
        Real& s2, Vec3& cp)
    {
        const Vec3 dq  = q1 - q0;
        const Real dqn = dq.norm2();

        if (dqn < s_eps)
        {
            s2 = Real(0);
            cp = q0;
            return true;
        }

        s2 = std::clamp((node - q0) * dq / dqn, Real(0), Real(1));
        cp = q0 + dq * s2;
        return true;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  candidateSegments – broad-phase bounding-sphere overlap test.
    //
    // These are abstract broad-phase radii supplied
    //  by the caller (contact-relevant: bore for outer tube in nested mode,
    //  outer radius otherwise). The implementation is unchanged.
    // ─────────────────────────────────────────────────────────────────────────────
    sofa::type::vector<int> SphereSweptIntersectionMethod::candidateSegments(
        const Vec3&          queryMid,
        Real                 queryHalfLength,
        Real                 r_bp_query,      
        const VecRigidCoord& candidateFrames,  
        Real                 r_bp_candidate,  
        Real                 broadPhaseMarginFactor)
    {
        const int N = static_cast<int>(candidateFrames.size());
        sofa::type::vector<int> candidates;
        candidates.reserve(std::max(0, N - 1));

        const Real R_i = queryHalfLength + r_bp_query * broadPhaseMarginFactor;

        for (int j = 0; j < N - 1; ++j)
        {
            const Vec3 q0    = candidateFrames[j].getCenter();
            const Vec3 q1    = candidateFrames[j + 1].getCenter();
            const Vec3 mid_j = (q0 + q1) * Real(0.5);
            const Real R_j   = (q1 - q0).norm() * Real(0.5)
                               + r_bp_candidate * broadPhaseMarginFactor;

            if ((queryMid - mid_j).norm() <= R_i + R_j)
                candidates.push_back(j);
        }
        return candidates;
    }

	// ─────────────────────────────────────────────────────────────────────────
    //  computeContactFrame
    //
    //  Algorithm:
    //    t̂₁ = normalize(τ₁ − (τ₁·n̂)·n̂)   [project onto contact plane]
    //    t̂₂ = n̂ × t̂₁
    //
    //  Fallback chain when τ₁ ∥ n̂:
    //    → try τ₂ (Beam-2 chord)
    //    → try global X, then global Y
    // ─────────────────────────────────────────────────────────────────────────
    void SphereSweptIntersectionMethod::computeContactFrame(
        const Vec3& tau1,
        const Vec3& tau2,
        const Vec3& nHat,
        Vec3&       t1_out,
        Vec3&       t2_out) 
    {
        // Project tau1 onto the contact plane (perpendicular to nHat).
        const Vec3 proj1  = tau1 - nHat * (tau1 * nHat);
        const Real proj1n = proj1.norm();
 
        if (proj1n > s_eps)
        {
            t1_out = proj1 / proj1n;
        }
        else
        {
            // tau1 ∥ nHat: try tau2.
            const Vec3 proj2  = tau2 - nHat * (tau2 * nHat);
            const Real proj2n = proj2.norm();
            if (proj2n > s_eps)
            {
                t1_out = proj2 / proj2n;
            }
            else
            {
                Vec3 ref     = Vec3(Real(1), Real(0), Real(0));
                Vec3 ref_proj = ref - nHat * (ref * nHat);
                if (ref_proj.norm() < s_eps)
                {
                    ref      = Vec3(Real(0), Real(1), Real(0));
                    ref_proj = ref - nHat * (ref * nHat);
                }
                const Real rn = ref_proj.norm();
                t1_out = (rn > s_eps) ? ref_proj / rn : Vec3(Real(1), Real(0), Real(0));
            }
        }
 
        // t̂₂ = n̂ × t̂₁  (circumferential, unit by construction when n̂ ⊥ t̂₁)
        t2_out = Vec3(
            nHat[1]*t1_out[2] - nHat[2]*t1_out[1],
            nHat[2]*t1_out[0] - nHat[0]*t1_out[2],
            nHat[0]*t1_out[1] - nHat[1]*t1_out[0]);
        const Real t2n = t2_out.norm();
        if (t2n > s_eps) t2_out /= t2n;
    }
    
    // ─────────────────────────────────────────────────────────────────────────  
    //  axialOverlapMaxRadial — nested-mode segment-pair contact metric.          
    //  See header for theory.                                                    
    // ─────────────────────────────────────────────────────────────────────────  
    bool SphereSweptIntersectionMethod::axialOverlapMaxRadial(                    
        const Vec3& p0, const Vec3& p1,                                           
        const Vec3& q0, const Vec3& q1,                                           
        Real& s1_out, Real& s2_out,                                               
        Vec3& cp1, Vec3& cp2,                                                     
        Real& radial)                                                             
    {                                                                             
        const Vec3 d1   = p1 - p0;                                                
        const Real d1n2 = d1.norm2();                                             
        if (d1n2 < s_eps) return false;  // Degenerate T1 segment.                
                                                                                  
        const Vec3 d2          = q1 - q0;                                         
        const Real d2_dot_d1   = d2 * d1;                                         
        const Vec3 q0_minus_p0 = q0 - p0;                                         
        const Real q0_dot_d1   = q0_minus_p0 * d1;                                
                                                                                  
        // s_1 parameter on T1 where q0 / q1 project axially.                     
        const Real s1_at_q0 = q0_dot_d1 / d1n2;                                   
        const Real s1_at_q1 = (q0_dot_d1 + d2_dot_d1) / d1n2;                     
                                                                                  
        // Axial overlap clamped to T1's parametric range [0,1].                  
        const Real s1_lo    = std::min(s1_at_q0, s1_at_q1);                       
        const Real s1_hi    = std::max(s1_at_q0, s1_at_q1);                       
        const Real s1_start = std::max(s1_lo, Real(0));                           
        const Real s1_end   = std::min(s1_hi, Real(1));                           
                                                                                  
        if (s1_start > s1_end) return false;  // No axial overlap.                
                                                                                  
        // Degeneracy: T2 perpendicular to T1 axis → s_2(s_1) undefined.          
        // Use T2 midpoint as a conservative single-point representative.         
        const Real d1n         = std::sqrt(d1n2);                                 
        const Real d2n         = d2.norm();                                       
        const bool d2_perp_d1  = (std::abs(d2_dot_d1) < s_eps * d1n * d2n);       
                                                                                  
        auto evalAtS1 = [&](Real s1, Real& s2, Vec3& P, Vec3& Q)                  
        {                                                                         
            P = p0 + d1 * s1;                                                     
            if (d2_perp_d1)                                                       
            {                                                                     
                s2 = Real(0.5);                                                   
                Q  = (q0 + q1) * Real(0.5);                                       
            }                                                                     
            else                                                                  
            {                                                                     
                s2 = (s1 * d1n2 - q0_dot_d1) / d2_dot_d1;                         
                s2 = std::clamp(s2, Real(0), Real(1));                            
                Q  = q0 + d2 * s2;                                                
            }                                                                     
        };                                                                        
                                                                                  
        // r²(s_1) convex quadratic on the overlap → max at one of the endpoints. 
        Real s2a, s2b;                                                            
        Vec3 Pa, Pb, Qa, Qb;                                                      
        evalAtS1(s1_start, s2a, Pa, Qa);                                          
        evalAtS1(s1_end,   s2b, Pb, Qb);                                          
                                                                                  
        const Real r_a = (Qa - Pa).norm();                                        
        const Real r_b = (Qb - Pb).norm();                                        
        
        const Real tol = std::max(Real(1e-12), Real(1e-9) * std::max(r_a, r_b));
        if (std::abs(r_a - r_b) < tol)
        {
            // Tie — collapse to the overlap midpoint to avoid step-to-step flipping.
            const Real s1_mid = Real(0.5) * (s1_start + s1_end);
            Real s2_mid; Vec3 Pm, Qm;
            evalAtS1(s1_mid, s2_mid, Pm, Qm);
            s1_out = s1_mid; s2_out = s2_mid;
            cp1 = Pm; cp2 = Qm;
            radial = (Qm - Pm).norm();
        }
        else if (r_a >= r_b)                                                           
        {                                                                         
            s1_out = s1_start; s2_out = s2a;                                      
            cp1 = Pa; cp2 = Qa;                                                   
            radial = r_a;                                                         
        }                                                                         
        else                                                                      
        {                                                                         
            s1_out = s1_end; s2_out = s2b;                                        
            cp1 = Pb; cp2 = Qb;                                                   
            radial = r_b;                                                         
        }                                                                         
        return true;                                                              
    }                                                                                                                                      
     
 

    SphereSweptIntersectionMethod::Vec3
    SphereSweptIntersectionMethod::getContactTangent1(std::size_t k) const
    {
        ensureUpdated();  
        if (k >= m_contactTangents1.size())
        {
            msg_error() << "getContactTangent1(): index " << k
                        << " out of range (size=" << m_contactTangents1.size()
                        << "). Returning fallback (1,0,0).";
            return Vec3(Real(1), Real(0), Real(0));
        }
        return m_contactTangents1[k];
    }
 
    SphereSweptIntersectionMethod::Vec3
    SphereSweptIntersectionMethod::getContactTangent2(std::size_t k) const
    {
        ensureUpdated();  
        if (k >= m_contactTangents2.size())
        {
            msg_error() << "getContactTangent2(): index " << k
                        << " out of range (size=" << m_contactTangents2.size()
                        << "). Returning fallback (0,1,0).";
            return Vec3(Real(0), Real(1), Real(0));
        }
        return m_contactTangents2[k];
    }
    
    SphereSweptIntersectionMethod::Real
    SphereSweptIntersectionMethod::gapSignForPublishedNormal() const
    {
        const bool isNested =
            (d_contactConfiguration.getValue().getSelectedItem() == "nested");
        if (!isNested)
            return Real(1);
        const bool beam1IsOuter = (d_radius1.getValue() >= d_radius2.getValue());
        return beam1IsOuter ? Real(-1) : Real(1);
    }
    
    std::size_t SphereSweptIntersectionMethod::getNumContacts() const
    {
        ensureUpdated();  
        return m_contactNormals.size();
    }

    SphereSweptIntersectionMethod::Vec2d
SphereSweptIntersectionMethod::getCurvilinearParams(std::size_t k) const
    {
        ensureUpdated();  
        const auto& params = d_curvilinearParams.getValue();
        if (k >= params.size())
        {
            msg_error() << "getCurvilinearParams: index " << k
                        << " out of range (size=" << params.size() << ").";
            return Vec2d(Real(0), Real(0));
        }
        return params[k];
    }

    SphereSweptIntersectionMethod::Vec3
    SphereSweptIntersectionMethod::getDistances(std::size_t k) const
    {
        ensureUpdated();  
        if (k >= m_distances.size())
        {
            msg_error() << "getDistances: index " << k
                        << " out of range (size=" << m_distances.size() << ").";
            return Vec3(Real(0), Real(0), Real(0));
        }
        return m_distances[k];
    }

    SphereSweptIntersectionMethod::Vec3
    SphereSweptIntersectionMethod::getCenterlinePoint1(std::size_t k) const
    {
        ensureUpdated();  
        if (k >= m_centerlinePoints1.size())
        {
            msg_error() << "getCenterlinePoint1: index " << k
                        << " out of range (size=" << m_centerlinePoints1.size() << ").";
            return Vec3(Real(0), Real(0), Real(0));
        }
        return m_centerlinePoints1[k];
    }

    SphereSweptIntersectionMethod::Vec3
    SphereSweptIntersectionMethod::getCenterlinePoint2(std::size_t k) const
    {
        ensureUpdated();  
        if (k >= m_centerlinePoints2.size())
        {
            msg_error() << "getCenterlinePoint2: index " << k
                        << " out of range (size=" << m_centerlinePoints2.size() << ").";
            return Vec3(Real(0), Real(0), Real(0));
        }
        return m_centerlinePoints2[k];
    }

    SphereSweptIntersectionMethod::Vec3
    SphereSweptIntersectionMethod::getSurfacePoint1(std::size_t k) const
    {
        ensureUpdated();  
        if (k >= m_surfacePoints1.size())
        {
            msg_error() << "getSurfacePoint1: index " << k
                        << " out of range (size=" << m_surfacePoints1.size() << ").";
            return Vec3(Real(0), Real(0), Real(0));
        }
        return m_surfacePoints1[k];
    }

    SphereSweptIntersectionMethod::Vec3
    SphereSweptIntersectionMethod::getSurfacePoint2(std::size_t k) const
    {   
        ensureUpdated();  
        if (k >= m_surfacePoints2.size())
        {
            msg_error() << "getSurfacePoint2: index " << k
                        << " out of range (size=" << m_surfacePoints2.size() << ").";
            return Vec3(Real(0), Real(0), Real(0));
        }
        return m_surfacePoints2[k];
    }

    SphereSweptIntersectionMethod::Vec2i
    SphereSweptIntersectionMethod::getContactSectionIds(std::size_t k) const
    {   
        ensureUpdated();  
        if (k >= m_contactSectionIds.size())
        {
            msg_error() << "getContactSectionIds: index " << k
                        << " out of range (size=" << m_contactSectionIds.size() << ").";
            return Vec2i(-1, -1);
        }
        return m_contactSectionIds[k];
    }
    
    void SphereSweptIntersectionMethod::ensureUpdated() const
    {
        // updateIfDirty() comes from DDGNode (via DataEngine). It is a no-op
        // if the node is already clean, otherwise it runs doUpdate().
        // const_cast is the documented SOFA idiom for this pattern: the logical
        // state (inputs → outputs) is unchanged, only the cached representation
        // is refreshed.
        const_cast<SphereSweptIntersectionMethod*>(this)->updateIfDirty();
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  SOFA factory registration
    // ─────────────────────────────────────────────────────────────────────────────
    void registerSphereSweptIntersectionMethod(sofa::core::ObjectFactory* factory)
    {
        factory->registerObjects(sofa::core::ObjectRegistrationData(
            "Computes the minimum distance and contact points between two "
            "Cosserat beams modelled as sphere-swept (canal) surfaces.\n"
            "Implements Lee 2007: minimum distance between two canal surfaces "
            "= minimum distance between two moving spheres along their centrelines.\n"
            "Supports solid beams (external contact) and hollow tubes (CTR nesting).")
            .add<SphereSweptIntersectionMethod>());
    }

} // namespace Cosserat
