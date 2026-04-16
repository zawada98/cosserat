/******************************************************************************
 * Cosserat Plugin for SOFA Framework                                         *
 *                                                                            *
 * BeamContactMapping.cpp                                                     *
 *                                                                            *
 * See BeamContactMapping.h for full documentation.                          *
 ******************************************************************************/
#include "BeamContactMapping.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/helper/accessor.h>
#include <sofa/type/Mat.h>
#include <algorithm>    // std::min

 // ── Explicit template instantiation ──────────────────────────────────────────
 // Multi2Mapping<TIn1,TIn2,TOut> is a class template whose method bodies live
 // in Multi2Mapping.inl.  No SOFA library ships a pre-built instantiation for
 // <Rigid3d, Rigid3d, Vec3d>, so the plugin must emit one here.
 //
 // Without this, the linker cannot resolve any base-class method:
 //   constructor, init(), apply(), applyJ(), applyJT() (both overloads),
 //   addInputModel1/2(), addOutputModel(), getFrom/To(), getMechFrom/To(),
 //   computeAccFromMapping(), disable().
 //
 // NOTE: do NOT add SOFA_CORE_API here.  Multi2Mapping originates from
 // Sofa.Core.dll; applying __declspec(dllexport) on an instantiation of a class
 // that was already __declspec(dllimport) in the same TU causes a
 // dllexport/dllimport conflict on MSVC.
#include <sofa/core/Multi2Mapping.inl>

namespace sofa::core
{
    template class Multi2Mapping<
        sofa::defaulttype::Rigid3dTypes,
        sofa::defaulttype::Rigid3dTypes,
        sofa::defaulttype::Vec3dTypes>;
}
// ─────────────────────────────────────────────────────────────────────────────

namespace Cosserat
{
    BeamContactMapping::BeamContactMapping()
        : Inherit1()
        , d_contactSectionIds(
            initData(&d_contactSectionIds,
                "contactSectionIds",
                "Per-contact {i,j} index pair from SphereSweptIntersectionMethod.\n"
                "ALGO_1: i = Beam-1 segment index, j = Beam-2 segment index.\n"
                "ALGO_2: i = Beam-1 node index,    j = Beam-2 segment index."))
        , d_curvilinearParams(
            initData(&d_curvilinearParams,
                "curvilinearParams",
                "Per-contact {alpha,beta} curvilinear parameters in [0,1] from SSIM.\n"
                "ALGO_2: alpha = 0 always (contact point coincides with Beam-1 node i)."))
        , d_radius1(
            initData(&d_radius1, Real(0.1),
                "radius1",
                "Cross-section radius of Beam-1 (same length unit as frame positions)."))
        , d_radius2(
            initData(&d_radius2, Real(0.1),
                "radius2",
                "Cross-section radius of Beam-2."))
        , d_isAlgo2(
            initData(&d_isAlgo2, false,
                "isAlgo2",
                "Set true when SSIM runs ALGO_2 (node-to-segment).\n"
                "Changes the role of sectionIds[k].x() from segment index to node index\n"
                "and forces nBeam1Blocks = 1 with weight = 1.")) //TODO: maybe change name
        , d_mappingMode(
            initData(&d_mappingMode, std::string("gap"),
                "mappingMode",
                // MODIFIED: both modes now require exactly ONE connected output MO.
                "Output mapping mode.  Both modes use exactly ONE connected output MO.\n"
                "  'contactPoints': output size = 2K (K = number of contact pairs).\n"
                "    Even indices : out[0][2k]   = Pc_A[k] = P_A + r1*n  (Beam-1 surface).\n"
                "    Odd  indices : out[0][2k+1] = Pc_B[k] = P_B - r2*n  (Beam-2 surface).\n"
                "    n = (P_B - P_A) / ||P_B - P_A||.\n"
                "    applyJ  gives interleaved velocities [Vc_A[0], Vc_B[0], Vc_A[1], ...].\n"
                "    applyJT back-projects: even inForce rows -> Beam-1, odd -> Beam-2.\n"
                "    applyJT(MatrixDeriv): even cols -> Beam-1, odd cols -> Beam-2.\n"
                "  'gap': output size = K.\n"
                "    out[0][k] = Pc_B[k] - Pc_A[k] = (d - r1 - r2)*n  (gap vector).\n"
                "    delta > 0: separation;  delta < 0: penetration.\n"
                "    applyJ  gives gap velocity: delta_dot = Pc_B_dot - Pc_A_dot.\n"
                "    applyJT back-projects: Beam-1 gets -w*F, Beam-2 gets +w*F."))
    {
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  init / reinit
    // ─────────────────────────────────────────────────────────────────────────────
    void BeamContactMapping::init()
    {
        // Validate the mode string at init time so the user gets a clear error
        // rather than silent misbehaviour at simulation start.
        std::cerr << ">>> [BCM] ctor done, nFields=" << this->getDataFields().size() << std::endl;

        const std::string& mode = d_mappingMode.getValue();
        if (mode != "contactPoints" && mode != "gap")
        {
            msg_error() << "Unknown mappingMode '" << mode << "'. "
                "Valid values are 'contactPoints' and 'gap'. "
                "Falling back to 'gap'.";
            d_mappingMode.setValue("gap");

        }

        Inherit1::init();
    }

    void BeamContactMapping::reinit()
    {
        Inherit1::reinit();
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  apply
    //
    //  Computes output positions and rebuilds the Jacobian cache.
    //
    //  ── Geometry (shared by both modes) ────────────────────────────────────────
    //
    //  ALGO_1 (isAlgo2 = false):
    //    P_A  = (1−α)·p[i]  + α·p[i+1]           (Beam-1 centreline)
    //    P_B  = (1−β)·p[j]  + β·p[j+1]           (Beam-2 centreline)
    //
    //  ALGO_2 (isAlgo2 = true), α = 0 always:
    //    P_A  = p[i]                               (Beam-1 node, no interpolation)
    //    P_B  = (1−β)·p[j]  + β·p[j+1]           (same as ALGO_1)
    //
    //  Common:
    //    n̂   = (P_B − P_A) / ‖P_B − P_A‖
    //    Pc_A = P_A + r₁·n̂
    //    Pc_B = P_B - r₂·n̂
    //
    //  ── output layout ──────────────────────────────────────────────────────────
    //
    //  mode = "contactPoints":
    //    Single output MO, size 2K (interleaved):
    //    out[0][2k]   = Pc_A[k]   (Beam-1 surface point)
    //    out[0][2k+1] = Pc_B[k]   (Beam-2 surface point)
    //
    //  mode = "gap":
    //    out[0][k] = Pc_B[k] − Pc_A[k]  =  (d − r₁ − r₂)·n̂
    //
    //  ── Jacobian cache ──────────────────────────────────────────────────────────
    //
    //  arm_A = +r₁·n̂  (contact-to-centreline vector for Beam-1 frames)
    //  arm_B = −r₂·n̂  (contact-to-centreline vector for Beam-2 frames)
    //
    //  Ṗc_A = Σ_b w_b · [v_b + ω_b × arm_A]
    //  Ṗc_B = Σ_b w_b · [v_b + ω_b × arm_B]
    //  δ̇    = Ṗc_B − Ṗc_A  (gap velocity, mode = "gap")
    // ─────────────────────────────────────────────────────────────────────────────
    void BeamContactMapping::apply(
        const sofa::core::MechanicalParams* /*mparams*/,
        const sofa::type::vector<sofa::core::objectmodel::Data<OutVecCoord>*>&
        dataVecOutPos,
        const sofa::type::vector<const sofa::core::objectmodel::Data<In1VecCoord>*>&
        dataVecIn1Pos,
        const sofa::type::vector<const sofa::core::objectmodel::Data<In2VecCoord>*>&
        dataVecIn2Pos)
    {
        if (dataVecIn1Pos.empty() || dataVecIn2Pos.empty() || dataVecOutPos.empty())
            return;

        sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<In1VecCoord>>
            frames1 = *dataVecIn1Pos[0];
        sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<In2VecCoord>>
            frames2 = *dataVecIn2Pos[0];

        const auto& ids = d_contactSectionIds.getValue();
        const auto& params = d_curvilinearParams.getValue();
        const Real  r1 = d_radius1.getValue();
        const Real  r2 = d_radius2.getValue();
        const bool  algo2 = d_isAlgo2.getValue();

        const int K = static_cast<int>(std::min(ids.size(), params.size()));

        if (ids.size() != params.size())
            msg_warning() << "contactSectionIds (" << ids.size()
            << ") and curvilinearParams (" << params.size()
            << ") have different sizes. Using K=" << K << ".";

        m_jacCache.resize(static_cast<sofa::Size>(K));

        const int N1 = static_cast<int>(frames1.size());
        const int N2 = static_cast<int>(frames2.size());

        // ── Per-contact geometry computation (shared by both modes) ──────────────
        //
        // Extracted into a lambda so it can be called identically from each branch
        // without duplicating the 60-line loop.  The lambda writes Pc_A/Pc_B and
        // the Jacobian cache entry; the caller decides which output slot to fill.
        //
        // Returns false if the contact pair is out of range (caller writes zeros).
        struct ContactGeom { Vec3 Pc_A, Pc_B; };

        auto computeGeom = [&](int k) -> std::pair<bool, ContactGeom>
            {
                const int  i = ids[k][0];
                const int  j = ids[k][1];
                const Real alpha = params[k][0];
                const Real beta = params[k][1];

                const bool b1Valid = algo2 ? (i >= 0 && i < N1)
                    : (i >= 0 && i + 1 < N1);
                const bool b2Valid = (j >= 0 && j + 1 < N2);

                if (!b1Valid || !b2Valid)
                {
                    msg_error() << "Contact pair " << k << ": index out of range "
                        << "(i=" << i << " j=" << j
                        << " N1=" << N1 << " N2=" << N2 << "). Skipping.";
                    m_jacCache[k].nBeam1Blocks = 0;
                    m_jacCache[k].nBeam2Blocks = 0;
                    return { false, {} };
                }

                // Beam-1 centreline point P_A
                Vec3 P_A;
                if (algo2)
                    P_A = frames1[i].getCenter();
                else
                    P_A = frames1[i].getCenter() * (Real(1) - alpha)
                    + frames1[i + 1].getCenter() * alpha;

                // Beam-2 centreline point P_B
                const Vec3 P_B = frames2[j].getCenter() * (Real(1) - beta)
                    + frames2[j + 1].getCenter() * beta;

                // Contact normal n̂ = (P_B − P_A) / ‖P_B − P_A‖
                Vec3 n = P_B - P_A;
                const Real dn = n.norm();
                if (dn < s_eps)
                {
                    Vec3 tangent = frames1[i].getOrientation().rotate(Vec3(1, 0, 0));
                    Vec3 cand[3] = { {1,0,0},{0,1,0},{0,0,1} }; //todo: change in later version
                    int bestIdx = 0;
                    for (int c = 1; c < 3; ++c)
                        if (std::abs(dot(cand[c], tangent)) < std::abs(dot(cand[bestIdx], tangent)))
                            bestIdx = c;
                    n = cross(tangent, cand[bestIdx]);
                    n.normalize();
                }
                else
                {
                    n /= dn;
                }

                // Surface points
                const Vec3 Pc_A = P_A + n * r1;   // Beam-1 surface toward Beam-2
                const Vec3 Pc_B = P_B - n * r2;   // Beam-2 surface toward Beam-1

                // Jacobian cache
                const Vec3 p_i = frames1[i].getCenter();
                const Vec3 p_j = frames2[j].getCenter();
                const Vec3 p_jp1 = frames2[j + 1].getCenter();

                const Vec3 armA1 = Pc_A - p_i;

                const Vec3 armB1 = Pc_B - p_j;
                const Vec3 armB2 = Pc_B - p_jp1;

                ContactJacEntry& entry = m_jacCache[k];
                if (algo2)
                {
                    entry.nBeam1Blocks = 1;
                    entry.beam1Blocks[0] = { i, Real(1), armA1 };
                }
                else
                {
                    const Vec3 p_ip1 = frames1[i + 1].getCenter();
                    const Vec3 armA2 = Pc_A - p_ip1;
                    entry.nBeam1Blocks = 2;
                    entry.beam1Blocks[0] = { i,     Real(1) - alpha, armA1 };
                    entry.beam1Blocks[1] = { i + 1, alpha,           armA2 };
                }
                entry.nBeam2Blocks = 2;
                entry.beam2Blocks[0] = { j,     Real(1) - beta, armB1 };
                entry.beam2Blocks[1] = { j + 1, beta,           armB2 };

                return { true, { Pc_A, Pc_B } };
            };

        // ── MONOTONIC RESIZE: never shrink the output MOs ────────────────────────
        //
        // storeLambda (called by FreeMotionAnimationLoop after the constraint solve)
        // writes back the resolved multiplier λ to each output MO using the DOF
        // indices that ContactFeeder registered at AnimateBeginEvent (indices 0..K_n-1).
        //
        // apply() is called a second time during the free-motion phase with the
        // integrated positions.  If SSIM detects K_free < K_n contact pairs in
        // the free-motion configuration, resizing to K_free would leave the MOs
        // smaller than what storeLambda expects: writing to index k ≥ K_free is
        // an out-of-bounds access → SIGSEGV.
        //
        // Fix: only grow the MOs, never shrink them.  Slots beyond the current K
        // retain stale data, but ContactFeeder only registers k < K_n contacts, so
        // storeLambda only accesses valid slots.  The ULC violation is only computed
        // for the active k range, so stale slots have no mechanical effect.

        // MODIFIED: contactPoints mode output size is now 2K (interleaved), not K.
        // Both modes share the single-output-MO path; only the slot count differs.
        const sofa::Size newK = static_cast<sofa::Size>(K);

        if (isGapMode())
        {
            // ── Gap mode: δ = Pc_B − Pc_A ──────────
            sofa::helper::WriteOnlyAccessor<sofa::core::objectmodel::Data<OutVecCoord>>
                outGap = *dataVecOutPos[0];
            if (outGap.size() < newK) outGap.resize(newK);

            for (int k = 0; k < K; ++k)
            {
                auto [ok, g] = computeGeom(k);
                if (!ok) continue;
                outGap[k] = g.Pc_B - g.Pc_A;
            }
        }
        else
        {
            // ── contactPoints mode: single interleaved output MO ─────────────────
            //    out[0][2k]   = Pc_A[k]   (even index → Beam-1 surface point)
            //    out[0][2k+1] = Pc_B[k]   (odd  index → Beam-2 surface point)
            sofa::helper::WriteOnlyAccessor<sofa::core::objectmodel::Data<OutVecCoord>>
                out = *dataVecOutPos[0];

            const sofa::Size newK2 = static_cast<sofa::Size>(2 * K);
            if (out.size() < newK2) out.resize(newK2);

            for (int k = 0; k < K; ++k)
            {
                auto [ok, g] = computeGeom(k);
                if (!ok) continue;
                out[2 * k]     = g.Pc_A;
                out[2 * k + 1] = g.Pc_B;
            }
        }
    }


    // ─────────────────────────────────────────────────────────────────────────────
    //  applyJ  –  velocity propagation
    //
    //  For each contributing frame block b:
    //    Ṗc[contactPoint] += w_b · (v_b + ω_b × arm_b)
    //
    //  mode = "contactPoints":
    //    Single output velocity vector, interleaved:
    //    outVel[0][2k]   = Ṗc_A[k]  (Beam-1 frames, arm = +r₁·n̂)
    //    outVel[0][2k+1] = Ṗc_B[k]  (Beam-2 frames, arm = −r₂·n̂)
    //
    //  mode = "gap":
    //    outVel[0][k] = δ̇[k] = Ṗc_B[k] − Ṗc_A[k]
    // ─────────────────────────────────────────────────────────────────────────────
    void BeamContactMapping::applyJ(
        const sofa::core::MechanicalParams* /*mparams*/,
        const sofa::type::vector<sofa::core::objectmodel::Data<OutVecDeriv>*>&
        dataVecOutVel,
        const sofa::type::vector<const sofa::core::objectmodel::Data<In1VecDeriv>*>&
        dataVecIn1Vel,
        const sofa::type::vector<const sofa::core::objectmodel::Data<In2VecDeriv>*>&
        dataVecIn2Vel)
    {
        if (dataVecIn1Vel.empty() || dataVecIn2Vel.empty() || dataVecOutVel.empty())
            return;

        sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<In1VecDeriv>>
            vel1 = *dataVecIn1Vel[0];
        sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<In2VecDeriv>>
            vel2 = *dataVecIn2Vel[0];

        const int  K = static_cast<int>(m_jacCache.size()); //todo: why not computed from calling params?
        const sofa::Size newK = static_cast<sofa::Size>(K);
        const bool gapMode = isGapMode();

        if (gapMode)
        {
            // ── Gap mode: single output, δ̇[k] = Ṗc_B[k] − Ṗc_A[k] ───────────────
            sofa::helper::WriteOnlyAccessor<sofa::core::objectmodel::Data<OutVecDeriv>>
                outVel = *dataVecOutVel[0];
            if (outVel.size() < newK) outVel.resize(newK);

            for (int k = 0; k < K; ++k)
            {
                const ContactJacEntry& entry = m_jacCache[k];

                OutDeriv vcA{};
                for (int b = 0; b < entry.nBeam1Blocks; ++b)
                {
                    const JacBlock& blk = entry.beam1Blocks[b];
                    vcA += (vel1[blk.frameIdx].getLinear()
                        + sofa::type::cross(vel1[blk.frameIdx].getAngular(),
                            blk.arm))
                        * blk.weight;
                }

                OutDeriv vcB{};
                for (int b = 0; b < entry.nBeam2Blocks; ++b)
                {
                    const JacBlock& blk = entry.beam2Blocks[b];
                    vcB += (vel2[blk.frameIdx].getLinear()
                        + sofa::type::cross(vel2[blk.frameIdx].getAngular(),
                            blk.arm))
                        * blk.weight;
                }

                outVel[k] = vcB - vcA;   // δ̇[k] = Ṗc_B − Ṗc_A
            }
        }
        else
        {
            // ── contactPoints mode: single interleaved output ─────────────────────
            //    outVel[0][2k]   = Ṗc_A[k]
            //    outVel[0][2k+1] = Ṗc_B[k]
            sofa::helper::WriteOnlyAccessor<sofa::core::objectmodel::Data<OutVecDeriv>>
                outVel = *dataVecOutVel[0];

            const sofa::Size newK2 = static_cast<sofa::Size>(2 * K);
            if (outVel.size() < newK2) outVel.resize(newK2);
            for (int k = 0; k < K; ++k)
            {
                const ContactJacEntry& entry = m_jacCache[k];

                OutDeriv vcA{};
                for (int b = 0; b < entry.nBeam1Blocks; ++b)
                {
                    const JacBlock& blk = entry.beam1Blocks[b];
                    vcA += (vel1[blk.frameIdx].getLinear()
                        + sofa::type::cross(vel1[blk.frameIdx].getAngular(),
                            blk.arm))
                        * blk.weight;
                }

                OutDeriv vcB{};
                for (int b = 0; b < entry.nBeam2Blocks; ++b)
                {
                    const JacBlock& blk = entry.beam2Blocks[b];
                    vcB += (vel2[blk.frameIdx].getLinear()
                        + sofa::type::cross(vel2[blk.frameIdx].getAngular(),
                            blk.arm))
                        * blk.weight;
                }

                outVel[2 * k]     = vcA;
                outVel[2 * k + 1] = vcB;
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  applyJT  (VecDeriv)  –  force back-propagation
    //
    //  Virtual-work principle:  δW = F · δPc = Jᵀ · F
    //
    //  mode = "contactPoints":                                                    // MODIFIED
    //    Single interleaved input force vector (size 2K):                         // MODIFIED
    //    inForce[0][2k]   = FA at Pc_A[k]  →  Beam-1 frames:                    // MODIFIED
    //      f += w_b·FA,   τ += w_b·(arm_A × FA)                                  // MODIFIED
    //    inForce[0][2k+1] = FB at Pc_B[k]  →  Beam-2 frames:                    // MODIFIED
    //      f += w_b·FB,   τ += w_b·(arm_B × FB)                                  // MODIFIED
    //
    //  mode = "gap":
    //    inForce[0][k] = F at gap DOF  →  both beams:
    //      Beam-1:  f -= w_b·F,   τ -= w_b·(arm_A × F)   (Ṗc_A enters with − in δ̇)
    //      Beam-2:  f += w_b·F,   τ += w_b·(arm_B × F)
    // ─────────────────────────────────────────────────────────────────────────────
    void BeamContactMapping::applyJT(
        const sofa::core::MechanicalParams* /*mparams*/,
        const sofa::type::vector<sofa::core::objectmodel::Data<In1VecDeriv>*>&
        dataVecOut1Force,
        const sofa::type::vector<sofa::core::objectmodel::Data<In2VecDeriv>*>&
        dataVecOut2Force,
        const sofa::type::vector<const sofa::core::objectmodel::Data<OutVecDeriv>*>&
        dataVecInForce)
    {
        if (dataVecInForce.empty() || dataVecOut1Force.empty() || dataVecOut2Force.empty())
            return;

        sofa::helper::WriteAccessor<sofa::core::objectmodel::Data<In1VecDeriv>>
            outF1 = *dataVecOut1Force[0];
        sofa::helper::WriteAccessor<sofa::core::objectmodel::Data<In2VecDeriv>>
            outF2 = *dataVecOut2Force[0];

        const int  K = static_cast<int>(m_jacCache.size());
        const bool gapMode = isGapMode();

        if (gapMode)
        {
            // Gap mode: single force input; back-project to both beams.
            sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<OutVecDeriv>>
                inForce = *dataVecInForce[0];

            for (int k = 0; k < K; ++k)
            {
                if (k >= static_cast<int>(inForce.size())) break;

                const ContactJacEntry& entry = m_jacCache[k];
                const Vec3 F = inForce[k];

                // Beam-1 receives negative contribution (Ṗc_A has − sign in δ̇).
                for (int b = 0; b < entry.nBeam1Blocks; ++b)
                {
                    const JacBlock& blk = entry.beam1Blocks[b];
                    outF1[blk.frameIdx].getLinear()
                        -= F * blk.weight;
                    outF1[blk.frameIdx].getAngular()
                        -= sofa::type::cross(blk.arm, F) * blk.weight;
                }

                // Beam-2 receives positive contribution.
                for (int b = 0; b < entry.nBeam2Blocks; ++b)
                {
                    const JacBlock& blk = entry.beam2Blocks[b];
                    outF2[blk.frameIdx].getLinear()
                        += F * blk.weight;
                    outF2[blk.frameIdx].getAngular()
                        += sofa::type::cross(blk.arm, F) * blk.weight;
                }
            }
        }
        else
        {
            // contactPoints mode: single interleaved force input.
            // inForce[0][2k]   = FA at Pc_A[k]  → Beam-1 frames only.
            // inForce[0][2k+1] = FB at Pc_B[k]  → Beam-2 frames only.
            sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<OutVecDeriv>>
                inForce = *dataVecInForce[0];

            for (int k = 0; k < K; ++k)
            {
                const ContactJacEntry& entry = m_jacCache[k];

                // Force at Pc_A[k]: even slot → Beam-1.
                const int slotA = 2 * k;
                if (slotA < static_cast<int>(inForce.size()))
                {
                    const Vec3 FA = inForce[slotA];
                    for (int b = 0; b < entry.nBeam1Blocks; ++b)
                    {
                        const JacBlock& blk = entry.beam1Blocks[b];
                        outF1[blk.frameIdx].getLinear()
                            += FA * blk.weight;
                        outF1[blk.frameIdx].getAngular()
                            += sofa::type::cross(blk.arm, FA) * blk.weight;
                    }
                }

                // Force at Pc_B[k]: odd slot → Beam-2.
                const int slotB = 2 * k + 1;
                if (slotB < static_cast<int>(inForce.size()))
                {
                    const Vec3 FB = inForce[slotB];
                    for (int b = 0; b < entry.nBeam2Blocks; ++b)
                    {
                        const JacBlock& blk = entry.beam2Blocks[b];
                        outF2[blk.frameIdx].getLinear()
                            += FB * blk.weight;
                        outF2[blk.frameIdx].getAngular()
                            += sofa::type::cross(blk.arm, FB) * blk.weight;
                    }
                }
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  applyJT  (MatrixDeriv)  –  constraint Jacobian assembly
    //
    //  Called by GenericConstraintSolver during the free-motion phase to build
    //  the per-constraint Jacobian rows that will be passed to the LCP/QP solver.
    //
    //  mode = "contactPoints":
    //    Single input constraint matrix (from the one output MO).
    //    Column index encoding:  col = 2k   → Pc_A[k] constraint → Beam-1 only.
    //                            col = 2k+1 → Pc_B[k] constraint → Beam-2 only.
    //    In1 frame b (Beam-1): w_b·d (translational), w_b·(arm_A × d) (rotational).
    //    In2 frame b (Beam-2): w_b·d (translational), w_b·(arm_B × d) (rotational).
    //
    //  mode = "gap":
    //    dataMatIn[0] holds gap constraint rows   → contributed to BOTH matrices.
    //      col k = contact index, d = constraint direction at gap DOF k.
    //      In1 frame b: −w_b·d (translational), −w_b·(arm_A × d) (rotational).
    //      In2 frame b: +w_b·d (translational), +w_b·(arm_B × d) (rotational).
    // ─────────────────────────────────────────────────────────────────────────────
    void BeamContactMapping::applyJT(
        const sofa::core::ConstraintParams* /*cparams*/,
        const sofa::type::vector<sofa::core::objectmodel::Data<In1MatrixDeriv>*>&
        dataMatOut1,
        const sofa::type::vector<sofa::core::objectmodel::Data<In2MatrixDeriv>*>&
        dataMatOut2,
        const sofa::type::vector<const sofa::core::objectmodel::Data<OutMatrixDeriv>*>&
        dataMatIn)
    {
        if (dataMatIn.empty() || dataMatOut1.empty() || dataMatOut2.empty())
            return;

        sofa::helper::WriteAccessor<sofa::core::objectmodel::Data<In1MatrixDeriv>>
            outM1 = *dataMatOut1[0];
        sofa::helper::WriteAccessor<sofa::core::objectmodel::Data<In2MatrixDeriv>>
            outM2 = *dataMatOut2[0];

        const int  K = static_cast<int>(m_jacCache.size());
        const bool gapMode = isGapMode();

        if (gapMode)
        {
            // Gap mode: single input constraint matrix.
            // Each column k = gap DOF index; contributions go to BOTH In1 and In2.
            sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<OutMatrixDeriv>>
                inMat = *dataMatIn[0];

            for (auto rowIt = inMat->begin(); rowIt != inMat->end(); ++rowIt)
            {
                typename In1MatrixDeriv::RowIterator row1 =
                    outM1->writeLine(rowIt.index());
                typename In2MatrixDeriv::RowIterator row2 =
                    outM2->writeLine(rowIt.index());

                for (auto colIt = rowIt.begin(); colIt != rowIt.end(); ++colIt)
                {
                    const int k = static_cast<int>(colIt.index());
                    if (k < 0 || k >= K) continue;

                    const Vec3 d = colIt.val();
                    const ContactJacEntry& entry = m_jacCache[k];

                    // Beam-1: negative Jacobian block (Ṗc_A enters with − in δ̇).
                    for (int b = 0; b < entry.nBeam1Blocks; ++b)
                    {
                        const JacBlock& blk = entry.beam1Blocks[b];
                        In1Deriv contrib;
                        contrib.getLinear() = -d * blk.weight;
                        contrib.getAngular() = -sofa::type::cross(blk.arm, d) * blk.weight;
                        row1.addCol(blk.frameIdx, contrib);
                    }

                    // Beam-2: positive Jacobian block.
                    for (int b = 0; b < entry.nBeam2Blocks; ++b)
                    {
                        const JacBlock& blk = entry.beam2Blocks[b];
                        In2Deriv contrib;
                        contrib.getLinear() = d * blk.weight;
                        contrib.getAngular() = sofa::type::cross(blk.arm, d) * blk.weight;
                        row2.addCol(blk.frameIdx, contrib);
                    }
                }
            }
        }
        else
        {
            // contactPoints mode: single interleaved constraint matrix.
            //   Even column 2k   → Pc_A[k] constraint → Beam-1 only.
            //   Odd  column 2k+1 → Pc_B[k] constraint → Beam-2 only.
            sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<OutMatrixDeriv>>
                inMat = *dataMatIn[0];

            for (auto rowIt = inMat->begin(); rowIt != inMat->end(); ++rowIt)
            {
                typename In1MatrixDeriv::RowIterator row1 =
                    outM1->writeLine(rowIt.index());
                typename In2MatrixDeriv::RowIterator row2 =
                    outM2->writeLine(rowIt.index());

                for (auto colIt = rowIt.begin(); colIt != rowIt.end(); ++colIt)
                {
                    const int col = static_cast<int>(colIt.index());
                    const int k   = col / 2;
                    if (k < 0 || k >= K) continue;

                    const Vec3 d = colIt.val();
                    const ContactJacEntry& entry = m_jacCache[k];

                    if (col % 2 == 0)
                    {
                        // Even column → Pc_A[k] → Beam-1 frames only.
                        for (int b = 0; b < entry.nBeam1Blocks; ++b)
                        {
                            const JacBlock& blk = entry.beam1Blocks[b];
                            In1Deriv contrib;
                            contrib.getLinear()  = d * blk.weight;
                            contrib.getAngular() = sofa::type::cross(blk.arm, d) * blk.weight;
                            row1.addCol(blk.frameIdx, contrib);
                        }
                    }
                    else
                    {
                        // Odd column → Pc_B[k] → Beam-2 frames only.
                        for (int b = 0; b < entry.nBeam2Blocks; ++b)
                        {
                            const JacBlock& blk = entry.beam2Blocks[b];
                            In2Deriv contrib;
                            contrib.getLinear()  = d * blk.weight;
                            contrib.getAngular() = sofa::type::cross(blk.arm, d) * blk.weight;
                            row2.addCol(blk.frameIdx, contrib);
                        }
                    }
                }
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  SOFA factory registration
    // ─────────────────────────────────────────────────────────────────────────────
    void registerBeamContactMapping(sofa::core::ObjectFactory* factory)
    {
        std::cerr << ">>> [BCM] registerBeamContactMapping called" << std::endl;
        factory->registerObjects(sofa::core::ObjectRegistrationData(
            "Maps SSIM contact-point descriptors (section IDs + curvilinear parameters) "
            "to two sets of Cosserat beam Rigid3d frames.\n"
            // MODIFIED: description updated for single-MO contactPoints layout.
            "Both modes use exactly ONE connected output MechanicalObject.\n"
            "Selectable via mappingMode:\n"
            "  'contactPoints': output size = 2K (interleaved). "
            "out[0][2k]=Pc_A (Beam-1 surface), out[0][2k+1]=Pc_B (Beam-2 surface). "
            "applyJT: even force/column indices -> Beam-1, odd -> Beam-2.\n"
            "  'gap': output size = K. "
            "out[0][k] = Pc_B[k]-Pc_A[k] = (d-r1-r2)*n (gap vector). "
            "applyJ gives gap velocity; applyJT back-projects with opposite signs to each beam.\n"
            "Supports ALGO_1 (segment-to-segment) and ALGO_2 (node-to-segment).\n"
            "Implements apply / applyJ / applyJT(VecDeriv) / applyJT(MatrixDeriv) "
            "for use with FreeMotionAnimationLoop + GenericConstraintSolver.")
            .add<BeamContactMapping>());
    }

} // namespace Cosserat