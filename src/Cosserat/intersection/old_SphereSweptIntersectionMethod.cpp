/******************************************************************************
 * Cosserat Plugin for SOFA Framework                                         *
 *                                                                            *
 * SphereSweptIntersectionMethod.cpp                                          *
 *                                                                            *
 * See SphereSweptIntersectionMethod.h for full documentation.               *
 ******************************************************************************/
#include "SphereSweptIntersectionMethod.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/accessor.h>
#include <sofa/type/Mat.h>          // Mat<3,3,Real> → needed by fromMatrix()

namespace Cosserat {

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
        , d_radius1(initData(&d_radius1, Real(0.1),
            "radius1",
            "Cross-section radius of Beam 1."))
        , d_radius2(initData(&d_radius2, Real(0.1),
            "radius2",
            "Cross-section radius of Beam 2."))
        , d_algorithmType(initData(&d_algorithmType,
            sofa::helper::OptionsGroup({ "ALGO_1", "ALGO_2" }),
            "algorithmType",
            "Contact detection algorithm. ALGO_1 = segment-to-segment (Lee 2007). "
            "ALGO_2 = node-to-segment (Xun, N+1 queries per step)."))
        , d_maxNRIterations(initData(&d_maxNRIterations, int(20),
            "maxNRIterations",
            "Maximum Newton-Raphson iterations (ALGO_2 only)."))
        , d_nrTolerance(initData(&d_nrTolerance, Real(1e-12),
            "nrTolerance",
            "Newton-Raphson convergence tolerance (ALGO_2 only)."))
        //
        // ── Outputs ───────────────────────────────────────────────────────────
        , d_curvilinearParams(initData(&d_curvilinearParams,
            "curvilinearParams",
            "Normalised curvilinear parameters {s1*, s2*} for each contact pair. "
            "ALGO_2: s1 = 0 always (contact at Beam-1 node)."))
        , d_distances(initData(&d_distances,
            "distances",
            "Gap vector {delta_n, delta_t1, delta_t2}. "
            "delta_n = ||Pint1-Pint2|| - (r1+r2), negative => penetration. "
            "Tangential components are reserved and set to zero."))
        , d_centerlinePoints1(initData(&d_centerlinePoints1,
            "centerlinePoints1",
            "Closest point on Beam 1 centreline for each contact pair."))
        , d_centerlinePoints2(initData(&d_centerlinePoints2,
            "centerlinePoints2",
            "Closest point on Beam 2 centreline for each contact pair."))
        , d_surfacePoints1(initData(&d_surfacePoints1,
            "surfacePoints1",
            "Surface contact point on Beam 1 as Vec3d. "
            "Position = Pint1 + r1*nhat  (nhat points Beam-1 -> Beam-2)."))
        , d_surfacePoints2(initData(&d_surfacePoints2,
            "surfacePoints2",
            "Surface contact point on Beam 2 as Vec3d. "
            "Position = Pint2 - r2*nhat  (nhat points Beam-1 -> Beam-2)."))
        , d_contactSectionIds(initData(&d_contactSectionIds,
            "contactSectionIds",
            "Beam-section (ALGO_1) or node/segment (ALGO_2) index pair {i, j} "
            "for each contact."))
    {
        sofa::helper::OptionsGroup algoOptions({ "ALGO_1", "ALGO_2" });
        algoOptions.setSelectedItem(0u); // default: ALGO_1
        d_algorithmType.setValue(algoOptions);
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
        addInput(&d_algorithmType);
        addInput(&d_maxNRIterations);
        addInput(&d_nrTolerance);

        addOutput(&d_curvilinearParams);
        addOutput(&d_distances);
        addOutput(&d_centerlinePoints1);
        addOutput(&d_centerlinePoints2);
        addOutput(&d_surfacePoints1);
        addOutput(&d_surfacePoints2);
        addOutput(&d_contactSectionIds);

        setDirtyValue();
    }

    void SphereSweptIntersectionMethod::reinit()
    {
        update();
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  doUpdate  – main entry point called every time inputs are dirty
    // ─────────────────────────────────────────────────────────────────────────────
    void SphereSweptIntersectionMethod::doUpdate()
    {
        // ── Read inputs ──────────────────────────────────────────────────────────
        sofa::helper::ReadAccessor<Data<VecRigidCoord>> frames1 = d_beam1Frames;
        sofa::helper::ReadAccessor<Data<VecRigidCoord>> frames2 = d_beam2Frames;

        const Real r1 = d_radius1.getValue();
        const Real r2 = d_radius2.getValue();
        const int  maxIter = d_maxNRIterations.getValue();
        const Real tol = d_nrTolerance.getValue();
        const bool useAlgo2 =
            (d_algorithmType.getValue().getSelectedItem() == std::string("ALGO_2"));

        const int N1 = static_cast<int>(frames1.size());
        const int N2 = static_cast<int>(frames2.size());

        if (N1 < 2 || N2 < 2)
        {
            msg_warning() << "Need at least 2 frames per beam. "
                "Beam1 has " << N1 << ", Beam2 has " << N2 << ".";
            return;
        }

        // ── Write-accessors for outputs ──────────────────────────────────────────
        sofa::helper::WriteOnlyAccessor<Data<sofa::type::vector<sofa::type::Vec2d>>>
            outParams = d_curvilinearParams;
        sofa::helper::WriteOnlyAccessor<Data<VecVec3>>       outDist = d_distances;
        sofa::helper::WriteOnlyAccessor<Data<VecVec3>>       outCL1 = d_centerlinePoints1;
        sofa::helper::WriteOnlyAccessor<Data<VecVec3>>       outCL2 = d_centerlinePoints2;
        sofa::helper::WriteOnlyAccessor<Data<VecVec3>> outS1 = d_surfacePoints1;
        sofa::helper::WriteOnlyAccessor<Data<VecVec3>> outS2 = d_surfacePoints2;
        sofa::helper::WriteOnlyAccessor<Data<VecVec2i>>      outIds = d_contactSectionIds;

        outParams.clear();
        outDist.clear();
        outCL1.clear();
        outCL2.clear();
        outS1.clear();
        outS2.clear();
        outIds.clear();

        // ─────────────────────────────────────────────────────────────────────────
        //  ALGO_1 – Segment-to-Segment
        //  Outer loop: N1-1 segments on Beam 1.
        //  Inner loop: all candidate segments on Beam 2 (broad-phase placeholder).
        // ─────────────────────────────────────────────────────────────────────────
        if (!useAlgo2)
        {
            for (int i = 0; i < N1 - 1; ++i)
            {
                const Vec3 p0 = frames1[i].getCenter();
                const Vec3 p1 = frames1[i + 1].getCenter();

                const auto candidates = candidateSections(i, frames1.ref(), frames2.ref());

                Real bestDist = std::numeric_limits<Real>::max();
                Real best_s1 = Real(0);
                Real best_s2 = Real(0);
                Vec3 best_cp1, best_cp2;
                int  best_j = -1;
                bool foundAny = false;

                for (const int j : candidates)
                {
                    const Vec3 q0 = frames2[j].getCenter();
                    const Vec3 q1 = frames2[j + 1].getCenter();

                    Real s1 = Real(0), s2 = Real(0);
                    Vec3 cp1, cp2;

                    if (!segmentToSegment(p0, p1, q0, q1, s1, s2, cp1, cp2))
                        continue;

                    const Real centrelineDist = (cp1 - cp2).norm();
                    const Real gap = centrelineDist - (r1 + r2);

                    if (gap < bestDist)
                    {
                        bestDist = gap;
                        best_s1 = s1;
                        best_s2 = s2;
                        best_cp1 = cp1;
                        best_cp2 = cp2;
                        best_j = j;
                        foundAny = true;
                    }
                }

                if (!foundAny) continue;

                Vec3 psurf1, psurf2;
               
                mapToSurface(best_cp1, best_cp2,                        
                    frames1[i],
                    r1, r2,
                    psurf1, psurf2);

                outParams.push_back({ best_s1, best_s2 });
                outDist.push_back(Vec3(bestDist, Real(0), Real(0))); //todo: add tangential components
                outCL1.push_back(best_cp1);
                outCL2.push_back(best_cp2);
                outS1.push_back(psurf1);
                outS2.push_back(psurf2);
                outIds.push_back({ i, best_j });
            }
        }
        // ─────────────────────────────────────────────────────────────────────────
        //  ALGO_2 – Node-to-Segment
        //  outer loop iterates over the N1 frame NODES on Beam 1,
        //  not over N1-1 sections with 3 probes each.
        //   This gives exactly N1 contact queries (= N_sections + 1).   
        //
        //  For each Beam-1 node i:
        //    s1 = 0  by definition (contact is located exactly at the node).
        //    Find the closest point on ALL Beam-2 segments → pick best j.
        // ─────────────────────────────────────────────────────────────────────────
        else
        {
            for (int i = 0; i < N1; ++i)
            {
                const Vec3 nodeP = frames1[i].getCenter();

                Real bestDist = std::numeric_limits<Real>::max();
                Real best_s2 = Real(0);
                Vec3 best_cp2;
                int  best_j = -1;
                bool foundAny = false;

                // Full scan over all Beam-2 segments (broad-phase placeholder)
                for (int j = 0; j < N2 - 1; ++j)
                {
                    const Vec3 q0 = frames2[j].getCenter();
                    const Vec3 q1 = frames2[j + 1].getCenter();

                    Real s2_nr = Real(0);
                    Vec3 cp2_nr;
                    if (!nodeToSegmentNR(nodeP, q0, q1, s2_nr, cp2_nr, maxIter, tol))
                        continue;

                    const Real centrelineDist = (nodeP - cp2_nr).norm();
                    const Real gap = centrelineDist - (r1 + r2);

                    if (gap < bestDist)
                    {
                        bestDist = gap;
                        best_s2 = s2_nr;
                        best_cp2 = cp2_nr;
                        best_j = j;
                        foundAny = true;
                    }
                }

                if (!foundAny) continue;

                Vec3 psurf1, psurf2;
                mapToSurface(nodeP, best_cp2,                            
                    frames1[i],
                    r1, r2,
                    psurf1, psurf2);

                outParams.push_back({ Real(0), best_s2 });
                outDist.push_back(Vec3(bestDist, Real(0), Real(0)));
                outCL1.push_back(nodeP);
                outCL2.push_back(best_cp2);
                outS1.push_back(psurf1);
                outS2.push_back(psurf2);
                outIds.push_back({ i, best_j });
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  mapToSurface
    //
    //  Computes surface contact positions (Vec3d) from centreline contact points.
    //  Previously: built two Rigid3d contact frames using SLERP-interpolated beam
    //  orientations and assembled a 3×3 rotation matrix per contact point.
    //  Now: position-only output; SLERP, Mat3, Quatd, and buildContactQuat removed.
    //
    //  psurf1 = pint1 + r1 * n̂   (on Beam 1 surface, pointing toward Beam 2)
    //  psurf2 = pint2 − r2 * n̂   (on Beam 2 surface, pointing toward Beam 1)
    //   where n̂ = (pint2 − pint1) / ||pint2 − pint1||   (Beam-1 → Beam-2)
    //
    //  Degenerate fallback when ||pint1−pint2|| < ε:
    //    n̂ = frameA local Y-axis (coincident centrelines).
    // ─────────────────────────────────────────────────────────────────────────────
    void SphereSweptIntersectionMethod::mapToSurface(
        const Vec3& pint1,
        const Vec3& pint2,
        const RigidCoord& frameA,
        Real r1, Real r2,
        Vec3& psurf1,
        Vec3& psurf2)
    {
        // ── 1. Contact normal n̂  (Pint1 → Pint2, Beam-1 → Beam-2) ─────────────
        Vec3 n = pint2 - pint1;  // modified: n̂ now points Beam-1 → Beam-2 (SOFA convention)
        const Real dnorm = n.norm();

        if (dnorm < s_eps)
        {
            // Degenerate: coincident centrelines → fall back to frameA local Y-axis
            n = frameA.getOrientation().rotate(Vec3(Real(0), Real(0),Real(1)));
            msg_warning("SphereSweptIntersectionMethod")
                << "Coincident centrelines detected. "
                "Falling back to frameA local z-axis as contact normal.";
        }
        else
        {
            n /= dnorm;
        }

        // ── 2. Surface positions ─────────────────────────────────────────────────
        psurf1 = pint1 + n * r1;  
        psurf2 = pint2 - n * r2;  
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  ALGO_1: Segment-to-Segment (Ericson 2005, §5.1.9)
    // ─────────────────────────────────────────────────────────────────────────────
    bool SphereSweptIntersectionMethod::segmentToSegment(
        const Vec3& p0, const Vec3& p1,
        const Vec3& q0, const Vec3& q1,
        Real& s1, Real& s2,
        Vec3& cp1, Vec3& cp2)
    {
        const Vec3 d1 = p1 - p0;   // direction of segment 1
        const Vec3 d2 = q1 - q0;   // direction of segment 2
        const Vec3 r = p0 - q0;

        const Real a = d1.norm2();  // squared length of seg 1
        const Real e = d2.norm2();  // squared length of seg 2
        const Real f = d2 * r;      // dot(d2, r)

        // ── Handle degenerate cases ──────────────────────────────────────────────
        if (a <= s_eps && e <= s_eps)
        {
            // Both segments degenerate to points
            s1 = s2 = Real(0);
            cp1 = p0;
            cp2 = q0;
            return true;
        }

        if (a <= s_eps)
        {
            // Segment 1 is a point
            s1 = Real(0);
            s2 = std::clamp(f / e, Real(0), Real(1));
        }
        else
        {
            const Real c = d1 * r;  // dot(d1, r)
            if (e <= s_eps)
            {
                // Segment 2 is a point
                s2 = Real(0);
                s1 = std::clamp(-c / a, Real(0), Real(1));
            }
            else
            {
                // General (non-degenerate) case
                const Real b = d1 * d2;          // dot(d1, d2)
                const Real denom = a * e - b * b;    // always >= 0

                if (denom > s_eps)
                    s1 = std::clamp((b * f - c * e) / denom, Real(0), Real(1));
                else
                    s1 = Real(0);  // parallel segments – pick one endpoint

                // Compute s2 from the unclamped result, then re-clamp
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
    //  ALGO_2: Node-to-Segment Newton–Raphson
    //
    //  Objective: min_s || P - Q(s) ||²  where Q(s) = q0 + s*(q1-q0)
    //
    //  Stationarity condition (Xun PhD, Eq. 3.8 / 3.9):
    //    f (s) = (Q(s) - P) · Q'(s) = 0
    //    f'(s) = Q'(s) · Q'(s)        (= const for linear segment)
    //
    //  NR update: s_{n+1} = s_n - f(s_n) / f'(s_n)
    //  Clamp s to [0,1] after each step.
    //
    // ─────────────────────────────────────────────────────────────────────────────
    bool SphereSweptIntersectionMethod::nodeToSegmentNR(
        const Vec3& node,
        const Vec3& q0, const Vec3& q1,
        Real& s2, Vec3& cp,
        int maxIter, Real tol)
    {
        const Vec3 dq = q1 - q0;          // Q'(s) – constant for a linear segment
        const Real dqn = dq.norm2();       // ||Q'||²

        if (dqn < s_eps)
        {
            // Degenerate segment (both endpoints coincide)
            s2 = Real(0);
            cp = q0;
            return true;
        }

        // Closed-form initial guess (exact for a linear segment in one iteration)
        s2 = std::clamp((node - q0) * dq / dqn, Real(0), Real(1));

        for (int iter = 0; iter < maxIter; ++iter)
        {
            const Vec3 Q_s = q0 + dq * s2;
            const Real f = (Q_s - node) * dq;   // stationarity residual

            if (std::abs(f) < tol * dqn)
                break;

            // NR step: f'(s) = dqn (constant)
            const Real ds = -f / dqn;
            const Real s2_new = std::clamp(s2 + ds, Real(0), Real(1));

            if (std::abs(s2_new - s2) < tol)
            {
                s2 = s2_new;
                break;
            }
            s2 = s2_new;
        }

        cp = q0 + dq * s2;
        return true;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  Broad-phase candidate selection
    //  currently returns all Beam-2 segments (no false negatives).
    // ─────────────────────────────────────────────────────────────────────────────
    sofa::type::vector<int>
        SphereSweptIntersectionMethod::candidateSections(
            int i,
            const VecRigidCoord& frames1,
            const VecRigidCoord& frames2)
    {
        const int N2 = static_cast<int>(frames2.size());
        sofa::type::vector<int> candidates;
        candidates.reserve(N2 - 1);

        const Vec3 mid1 = (frames1[i].getCenter() + frames1[i + 1].getCenter()) * Real(0.5);
        const Real segLen1 = (frames1[i + 1].getCenter() - frames1[i].getCenter()).norm();
        (void)mid1; (void)segLen1;   // suppress unused-variable warnings until BVH is added

        for (int j = 0; j < N2 - 1; ++j)
            candidates.push_back(j);

        return candidates;
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  SOFA factory registration
    // ─────────────────────────────────────────────────────────────────────────────
    void registerSphereSweptIntersectionMethod(sofa::core::ObjectFactory* factory)
    {
        factory->registerObjects(sofa::core::ObjectRegistrationData(
            "Computes the minimum distance and contact points between two "
            "Cosserat beams modelled as sphere-swept (canal) surfaces.\n"
            "Implements the Lee 2007 Lemma: the minimum distance between two "
            "canal surfaces equals the minimum distance between two moving "
            "spheres along their centrelines.")
            .add<SphereSweptIntersectionMethod>());
    }

} // namespace Cosserat