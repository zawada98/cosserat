/******************************************************************************
 * Cosserat Plugin for SOFA Framework                                         *
 *                                                                            *
 * BeamContactMapping.cpp                                                     *
 *                                                                            *
 * See BeamContactMapping.h for full documentation.                          *
 ******************************************************************************/
#include "BeamContactMapping.h"
#include "Cosserat/intersection/SphereSweptIntersectionMethod.h"  
#include <sofa/core/MechanicalParams.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/accessor.h>

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
#include <cstdlib>
#include <fstream>
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

    bool traceBeamContactMapping()
    {
        static const bool enabled = std::getenv("COSSERAT_TRACE_BCM") != nullptr;
        return enabled;
    }

    BeamContactMapping::BeamContactMapping()
        : Inherit1()
        , l_ssim(initLink("ssim",
                "Mandatory link to the SphereSweptIntersectionMethod that provides "
                "contact evaluations via evaluateContacts(). "
                "Set via the 'ssim' attribute in the SOFA scene, e.g.: "
                "ssim='@contact_node/ssim'."))
        , d_mappingMode(
            initData(&d_mappingMode, std::string("gap"),
                "mappingMode",
                "Output mapping mode.  Both modes use exactly ONE connected output MO.\n"
                "  'gap': output size = K.\n"
               "    out[0][k] = Vec3(delta_n, delta_t1, delta_t2) in contact-local frame.\n"
               "      delta_n  = (Pc_B - Pc_A).n  (signed normal gap, <0 = penetration).\n"
               "      delta_t1 = (Pc_B - Pc_A).t1 (axial tangential gap).\n"
               "      delta_t2 = (Pc_B - Pc_A).t2 (circumferential tangential gap).\n"
               "    Contact frame {n, t1, t2} matches SSIM ContactEvaluation distances:\n"
               "      t1 = normalize(tau1 - (tau1.n)*n), tau1 = Beam-1 segment chord.\n"
               "      t2 = n x t1.\n"
               "    applyJ  gives gap velocity Vec3(dPrel.n, dPrel.t1, dPrel.t2).\n"
               "    applyJT converts Vec3(F_n,F_t1,F_t2) -> F_phys = n*F_n + t1*F_t1 + t2*F_t2.\n"
               "  'contactPoints': output size = 2K (K = number of contact pairs).\n"
               "    Even indices : out[0][2k]   = Pc_A[k] = Pa + r1*n  (Beam-1 surface).\n"
               "    Odd  indices : out[0][2k+1] = Pc_B[k] = Pb - r2*n  (Beam-2 surface).\n"
               "    n = (Pb - Pa) / ||Pb - Pa||, fetched from SSIM.\n"
               "    applyJ  gives interleaved velocities [Vc_A[0], Vc_B[0], Vc_A[1], ...].\n"
               "    applyJT back-projects: even inForce rows -> Beam-1, odd -> Beam-2.\n"
               "    applyJT(MatrixDeriv): even cols -> Beam-1, odd cols -> Beam-2."))
        , d_contactTriads(
            initData(&d_contactTriads,        
                "contactTriads",
                "Per-pair contact triad (n̂, t̂₁, t̂₂) written by apply().\n"
                "  n   — unit contact normal, external → internal.\n"
                "  t1  — Beam-1 tangent projected onto the contact plane.\n"
                "  t2  — circumferential,  t2 = n × t1 (right-handed).\n"
                "Downstream constraints (CPULC) link to this field for both the\n"
                "normal row and, when μ > 0, the two friction rows."))
        , d_gapSign(initData(&d_gapSign, SReal(1),          
                "gapSign",
                "Global gap sign s ∈ {+1, −1} such that (Pc_B − Pc_A)·n̂ = s·δn.\n"
                "Fetched from SSIM::gapSignForPublishedNormal() in init().\n"
                "Downstream constraints link here to compute dfree consistently."))
        , d_distances(initData(&d_distances,
            "distances",
            "Consolidated gap Vec3(δn, δt1, δt2)[k] written by apply().\n"
            "δt1/δt2 come from SSIM ContactEvaluation distances[k][1,2] (velocity-integrated).\n"
            "Link DUCL/CPULC to this field instead of to SSIM."))
    {
    }

    // ─────────────────────────────────────────────────────────────────────────────
    //  init / reinit
    // ─────────────────────────────────────────────────────────────────────────────
    void BeamContactMapping::init()
    {
        // Validate SSIM link before anything else.
        if (!l_ssim.get())
        {
            msg_error() << "The 'ssim' link is not set. "
                           "BeamContactMapping requires a valid link to a "
                           "SphereSweptIntersectionMethod object to fetch contact normals. "
                           "Add ssim='@<path>/<ssimName>' to the addObject() call.";
            // Do not return: let the rest of init() run so SOFA reports all errors
            // at once rather than stopping at the first one.
        }

        const std::string& mode = d_mappingMode.getValue();
        if (mode != "contactPoints" && mode != "gap")
        {
            msg_error() << "Unknown mappingMode '" << mode << "'. "
                "Valid values are 'contactPoints' and 'gap'. "
                "Falling back to 'gap'.";
            d_mappingMode.setValue("gap");

        }
        
        if (l_ssim.get())
            d_gapSign.setValue(l_ssim->gapSignForPublishedNormal());
        
        Inherit1::init();
    }

    void BeamContactMapping::reinit()
    {
        Inherit1::reinit();
    }

    BeamContactMapping::EvaluationCacheKey BeamContactMapping::makeEvaluationCacheKey(
        const sofa::core::objectmodel::BaseData& frames1Data,
        const sofa::core::objectmodel::BaseData& frames2Data,
        const sofa::core::objectmodel::BaseData& vels1Data,
        const sofa::core::objectmodel::BaseData& vels2Data) const
    {
        EvaluationCacheKey key;
        key.ssim = l_ssim.get();
        key.ssimParameterCounter = l_ssim.get() ? l_ssim->getEvaluationParametersCounter() : 0;
        key.frames1 = &frames1Data;
        key.frames2 = &frames2Data;
        key.vels1 = &vels1Data;
        key.vels2 = &vels2Data;
        key.frames1Counter = frames1Data.getCounter();
        key.frames2Counter = frames2Data.getCounter();
        key.vels1Counter = vels1Data.getCounter();
        key.vels2Counter = vels2Data.getCounter();
        return key;
    }

    bool BeamContactMapping::isJacobianCacheValidFor(const EvaluationCacheKey& key) const
    {
        return m_jacCacheValid && key == m_jacCacheKey;
    }

    void BeamContactMapping::markJacobianCacheValidFor(const EvaluationCacheKey& key)
    {
        m_jacCacheKey = key;
        m_jacCacheValid = true;
    }

    bool BeamContactMapping::buildJacobianEntries(
        const SphereSweptIntersectionMethod::ContactEvaluation& eval,
        const In1VecCoord& frames1,
        const In2VecCoord& frames2,
        sofa::type::vector<ContactJacEntry>& entries,
        const char* caller) const
    {
        if (!l_ssim.get())
        {
            msg_error() << caller << "(): l_ssim is null. "
                           "Check that ssim='@...' is set in the scene.";
            entries.clear();
            return false;
        }

        const sofa::Size K = static_cast<sofa::Size>(eval.contactSectionIds.size());
        if (eval.curvilinearParams.size() < K ||
            eval.surfacePoints1.size() < K ||
            eval.surfacePoints2.size() < K ||
            eval.centerlinePoints1.size() < K ||
            eval.centerlinePoints2.size() < K ||
            eval.contactNormals.size() < K ||
            eval.contactTangents1.size() < K ||
            eval.contactTangents2.size() < K ||
            eval.distances.size() < K)
        {
            msg_error() << caller << "(): incomplete SSIM evaluation.";
            entries.clear();
            return false;
        }
        entries.resize(K);

        const int N1 = static_cast<int>(frames1.size());
        const int N2 = static_cast<int>(frames2.size());

        for (sofa::Size k = 0; k < K; ++k)
        {
            ContactJacEntry& entry = entries[k];

            const Vec2i sec = eval.contactSectionIds[k];
            const int   i   = sec[0];
            const int   j   = sec[1];
            const Vec2d cp  = eval.curvilinearParams[k];
            const Real  alpha = cp[0];
            const Real  beta  = cp[1];

            if (i < 0 || i + 1 >= N1 || j < 0 || j + 1 >= N2)
            {
                msg_error() << "Contact pair " << k << ": section index out of range "
                    << "(i=" << i << " j=" << j
                    << " N1=" << N1 << " N2=" << N2 << "). Skipping.";

                entry.normal   = Vec3(Real(0), Real(0), Real(0));
                entry.tangent1 = Vec3(Real(0), Real(0), Real(0));
                entry.tangent2 = Vec3(Real(0), Real(0), Real(0));
                entry.surfacePoint1 = Vec3(Real(0), Real(0), Real(0));
                entry.surfacePoint2 = Vec3(Real(0), Real(0), Real(0));
                entry.gapNormal   = s_invalidGap;
                entry.gapTangent1 = Real(0);
                entry.gapTangent2 = Real(0);
                continue;
            }

            const Vec3 Pc_A = eval.surfacePoints1[k];
            const Vec3 Pc_B = eval.surfacePoints2[k];
            const Vec3 P_A_ssim = eval.centerlinePoints1[k];
            const Vec3 P_B_ssim = eval.centerlinePoints2[k];

            entry.beam1Blocks[0] = { i,     Real(1) - alpha, Pc_A - P_A_ssim };
            entry.beam1Blocks[1] = { i + 1, alpha,           Pc_A - P_A_ssim };
            entry.beam2Blocks[0] = { j,     Real(1) - beta,  Pc_B - P_B_ssim };
            entry.beam2Blocks[1] = { j + 1, beta,            Pc_B - P_B_ssim };

            entry.surfacePoint1 = Pc_A;
            entry.surfacePoint2 = Pc_B;
            entry.normal        = eval.contactNormals[k];
            entry.tangent1      = eval.contactTangents1[k];
            entry.tangent2      = eval.contactTangents2[k];

            const Vec3 d = eval.distances[k];
            entry.gapNormal   = d[0];
            entry.gapTangent1 = d[1];
            entry.gapTangent2 = d[2];
        }

        return true;
    }

    bool BeamContactMapping::rebuildJacobianCache(
        const SphereSweptIntersectionMethod::ContactEvaluation& eval,
        const In1VecCoord& frames1,
        const In2VecCoord& frames2)
    {
        sofa::type::vector<ContactJacEntry> entries;
        if (!buildJacobianEntries(eval, frames1, frames2, entries, "rebuildJacobianCache"))
        {
            m_jacCache.clear();
            m_jacCacheValid = false;
            return false;
        }

        m_jacCache = entries;
        d_gapSign.setValue(l_ssim->gapSignForPublishedNormal());
        publishContactDataFromCache();
        return true;
    }

    bool BeamContactMapping::rebuildJacobianCacheForApplyJ(
        const sofa::core::MechanicalParams* mparams,
        const sofa::core::objectmodel::Data<In1VecDeriv>& vel1Data,
        const sofa::core::objectmodel::Data<In2VecDeriv>& vel2Data,
        sofa::type::vector<ContactJacEntry>& scratchCache,
        const sofa::type::vector<ContactJacEntry>*& jacCacheForApplyJ,
        Real& gapSignForApplyJ)
    {
        jacCacheForApplyJ = nullptr;
        scratchCache.clear();
        gapSignForApplyJ = d_gapSign.getValue();

        if (!l_ssim.get())
        {
            msg_error() << "rebuildJacobianCacheForApplyJ(): l_ssim is null.";
            return false;
        }

        const auto& from1 = this->getFromModels1();
        const auto& from2 = this->getFromModels2();
        if (from1.empty() || from2.empty() || !from1[0] || !from2[0])
        {
            msg_error() << "rebuildJacobianCacheForApplyJ(): missing input MechanicalState.";
            return false;
        }

        const auto* frames1Data = mparams
            ? from1[0]->read(mparams->x().getId(from1[0]))
            : nullptr;
        const auto* frames2Data = mparams
            ? from2[0]->read(mparams->x().getId(from2[0]))
            : nullptr;
        if (!frames1Data)
            frames1Data = from1[0]->read(sofa::core::vec_id::read_access::position);
        if (!frames2Data)
            frames2Data = from2[0]->read(sofa::core::vec_id::read_access::position);
        if (!frames1Data || !frames2Data)
        {
            msg_error() << "rebuildJacobianCacheForApplyJ(): could not read input positions.";
            return false;
        }

        const EvaluationCacheKey key =
            makeEvaluationCacheKey(*frames1Data, *frames2Data, vel1Data, vel2Data);
        if (isJacobianCacheValidFor(key))
        {
            jacCacheForApplyJ = &m_jacCache;
            gapSignForApplyJ = d_gapSign.getValue();
            return true;
        }

        const auto eval = l_ssim->evaluateContacts(
            *frames1Data, *frames2Data, vel1Data, vel2Data);
        if (!buildJacobianEntries(eval,
                                  frames1Data->getValue(),
                                  frames2Data->getValue(),
                                  scratchCache,
                                  "rebuildJacobianCacheForApplyJ"))
        {
            return false;
        }

        jacCacheForApplyJ = &scratchCache;
        gapSignForApplyJ = l_ssim->gapSignForPublishedNormal();
        return true;
    }

    bool BeamContactMapping::requireFrozenJacobianCache(const char* caller) const
    {
        if (m_jacCacheValid)
            return true;

        msg_error() << caller
                    << ": no frozen SSIM/BCM contact cache is available. "
                       "apply() must establish the contact cache before "
                       "transpose mapping is used. applyJ() intentionally "
                       "does not publish or overwrite the frozen cache.";
        return false;
    }

    void BeamContactMapping::publishContactDataFromCache()
    {
        auto triads = sofa::helper::getWriteOnlyAccessor(d_contactTriads);
        auto dists  = sofa::helper::getWriteOnlyAccessor(d_distances);
        const sofa::Size K = static_cast<sofa::Size>(m_jacCache.size());

        triads.resize(K);
        dists.resize(K);

        for (sofa::Size k = 0; k < K; ++k)
        {
            const ContactJacEntry& entry = m_jacCache[k];
            triads[k].n  = entry.normal;
            triads[k].t1 = entry.tangent1;
            triads[k].t2 = entry.tangent2;
            dists[k] = Vec3(entry.gapNormal, entry.gapTangent1, entry.gapTangent2);
        }
        
    }

    // ─────────────────────────────────────────────────────────────────────────
    //  apply
    //
    //  Computes output positions and rebuilds the Jacobian cache.
    //
    //  ── Geometry ─────────────────────────────────────────────────────────────
    //
    //  Surface contact points Pc_A and Pc_B are read from SSIM:       
    //    Pc_A = ContactEvaluation.surfacePoints1[k]
    //    Pc_B = ContactEvaluation.surfacePoints2[k]
    //  SSIM already accounts for all modes (external / nested CTR, solid / hollow)
    //  and applies the correct contact-relevant radii.  BCM never recomputes them.
    //
    //  Moment arms still use input-MO frame centres (Pc − p_frame), since
    //  those are not exported by SSIM.
    //
    //  ── Contact-local frame (gap mode) ───────────────────────────────────────  
    //
    //  Identical to SSIM ContactEvaluation distances convention:
    //    τ₁  = unit chord of Beam-1 segment [i, i+1]
    //    t̂₁ = normalize(τ₁ − (τ₁·n̂)·n̂)             (projected onto contact plane)
    //    t̂₂ = n̂ × t̂₁                                (circumferential)
    //
    //  Stored in m_jacCache[k].tangent1 / .tangent2 so applyJ / applyJT can
    //  use them without re-fetching from SSIM.
    //
    //  ── Output layout ────────────────────────────────────────────────────────
    //
    //  mode = "contactPoints":
    //    Single output MO, size 2K (interleaved):
    //    out[0][2k]   = Pc_A[k]   (Beam-1 surface point)
    //    out[0][2k+1] = Pc_B[k]   (Beam-2 surface point)
    //
    //  mode = "gap":                                                     
    //    out[0][k] = Vec3(δ_n, δ_t1, δ_t2)  (contact-local frame)
    //
    //  ── MONOTONIC RESIZE ─────────────────────────────────────────────────────
    //
    //  storeLambda writes back to the output MOs using DOF indices registered at
    //  AnimateBeginEvent.  If SSIM detects K_free < K_n during the free-motion
    //  phase and BCM shrinks the output MO, storeLambda accesses out-of-bounds
    //  slots → SIGSEGV.  Fix: only grow the MO, never shrink.  Stale slots beyond
    //  the active K have no mechanical effect.
    // ─────────────────────────────────────────────────────────────────────────
    void BeamContactMapping::apply(
        const sofa::core::MechanicalParams* mparams,
        const sofa::type::vector<sofa::core::objectmodel::Data<OutVecCoord>*>&
        dataVecOutPos,
        const sofa::type::vector<const sofa::core::objectmodel::Data<In1VecCoord>*>&
        dataVecIn1Pos,
        const sofa::type::vector<const sofa::core::objectmodel::Data<In2VecCoord>*>&
        dataVecIn2Pos)
    {
        if (traceBeamContactMapping())
        {
            static std::ofstream bcmlog("bcm_apply_log.txt", std::ios::out | std::ios::trunc);
            bcmlog << "[BCM.apply] t=" << this->getContext()->getTime()
                    << " xId=" << (mparams ? mparams->x().getName() : "null")
                   << " outPtr="  << static_cast<const void*>(dataVecOutPos[0])
                   << " in1Ptr="  << static_cast<const void*>(dataVecIn1Pos[0])
                   << " in2Ptr="  << static_cast<const void*>(dataVecIn2Pos[0])
                   << "\n";
            bcmlog.flush();
        }
        
        if (!l_ssim.get())
        {
            msg_error() << "apply(): l_ssim is null — cannot fetch SSIM outputs. "
                           "Aborting. Check that ssim='@...' is set in the scene.";
            return;
        }
        
        if (dataVecIn1Pos.empty() || dataVecIn2Pos.empty() || dataVecOutPos.empty())
            return;

        sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<In1VecCoord>>
            frames1 = *dataVecIn1Pos[0];
        sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<In2VecCoord>>
            frames2 = *dataVecIn2Pos[0];

        const auto& from1 = this->getFromModels1();
        const auto& from2 = this->getFromModels2();
        const auto* vels1Data = (!from1.empty() && from1[0] && mparams)
            ? from1[0]->read(mparams->v().getId(from1[0]))
            : nullptr;
        const auto* vels2Data = (!from2.empty() && from2[0] && mparams)
            ? from2[0]->read(mparams->v().getId(from2[0]))
            : nullptr; 
        if (!vels1Data && !from1.empty() && from1[0])
            vels1Data = from1[0]->read(sofa::core::vec_id::read_access::velocity);
        if (!vels2Data && !from2.empty() && from2[0])
            vels2Data = from2[0]->read(sofa::core::vec_id::read_access::velocity);
        if (!vels1Data || !vels2Data)
        {
            msg_error() << "apply(): could not read input velocities for SSIM evaluation.";
            return;
        }

        const EvaluationCacheKey key =
            makeEvaluationCacheKey(*dataVecIn1Pos[0], *dataVecIn2Pos[0],
                                   *vels1Data, *vels2Data);
        if (!isJacobianCacheValidFor(key))
        {
            const auto eval = l_ssim->evaluateContacts(
                *dataVecIn1Pos[0], *dataVecIn2Pos[0], *vels1Data, *vels2Data);
            if (!rebuildJacobianCache(eval, frames1.ref(), frames2.ref()))
                return;
            markJacobianCacheValidFor(key);
        }

        const sofa::Size cacheK = static_cast<sofa::Size>(m_jacCache.size());
        if (isGapMode())
        {
            sofa::helper::WriteOnlyAccessor<sofa::core::objectmodel::Data<OutVecCoord>>
                outGap = *dataVecOutPos[0];
            if (outGap.size() < cacheK) outGap.resize(cacheK);
            for (sofa::Size k = 0; k < cacheK; ++k)
            {
                const ContactJacEntry& entry = m_jacCache[k];
                outGap[k] = Vec3(entry.gapNormal,
                                 entry.gapTangent1,
                                 entry.gapTangent2);
            }
        }
        else
        {
            sofa::helper::WriteOnlyAccessor<sofa::core::objectmodel::Data<OutVecCoord>>
                out = *dataVecOutPos[0];
            const sofa::Size newK2 = 2 * cacheK;
            if (out.size() < newK2) out.resize(newK2);
            for (sofa::Size k = 0; k < cacheK; ++k)
            {
                const ContactJacEntry& entry = m_jacCache[k];
                out[2 * k]     = entry.surfacePoint1;
                out[2 * k + 1] = entry.surfacePoint2;
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
        const sofa::core::MechanicalParams* mparams,
        const sofa::type::vector<sofa::core::objectmodel::Data<OutVecDeriv>*>&
        dataVecOutVel,
        const sofa::type::vector<const sofa::core::objectmodel::Data<In1VecDeriv>*>&
        dataVecIn1Vel,
        const sofa::type::vector<const sofa::core::objectmodel::Data<In2VecDeriv>*>&
        dataVecIn2Vel)
    {
        if (dataVecIn1Vel.empty() || dataVecIn2Vel.empty() || dataVecOutVel.empty())
            return;

        sofa::type::vector<ContactJacEntry> scratchCache;
        const sofa::type::vector<ContactJacEntry>* jacCacheForApplyJ = nullptr;
        Real gapSignForApplyJ = Real(1);
        if (!rebuildJacobianCacheForApplyJ(mparams,
                                           *dataVecIn1Vel[0],
                                           *dataVecIn2Vel[0],
                                           scratchCache,
                                           jacCacheForApplyJ,
                                           gapSignForApplyJ) ||
            !jacCacheForApplyJ)
        {
            return;
        }

        if (traceBeamContactMapping())
        {
            static std::ofstream bcmlog("bcm_mapping_applyJ_log.txt",
                                       std::ios::out | std::ios::trunc);
            bcmlog << "[BCM.applyJ] t=" << this->getContext()->getTime()
                  << " mapping=" << this->getName()
                  << " dxId="      << (mparams ? mparams->dx().getName() : "null")
                  << " vId="       << (mparams ? mparams->v().getName()  : "null")
                  << " outPtr="  << static_cast<const void*>(dataVecOutVel[0])
                  << " in1Ptr="  << static_cast<const void*>(dataVecIn1Vel[0])
                  << " in2Ptr="  << static_cast<const void*>(dataVecIn2Vel[0])
                  << "\n";
            bcmlog.flush();
        }

        sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<In1VecDeriv>>
            vel1 = *dataVecIn1Vel[0];
        sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<In2VecDeriv>>
            vel2 = *dataVecIn2Vel[0];

        const auto& jacCache = *jacCacheForApplyJ;
        const sofa::Size  K = static_cast<int>(jacCache.size());
        const bool gapMode = isGapMode();
        const Real s = gapSignForApplyJ;
        
        if (gapMode)
        {
            // ── Gap mode: single output, δ̇[k] = Ṗc_B[k] − Ṗc_A[k] ───────────────
            sofa::helper::WriteOnlyAccessor<sofa::core::objectmodel::Data<OutVecDeriv>>
                outVel = *dataVecOutVel[0];
            if (outVel.size() < K) outVel.resize(K);

            for (sofa::Size k = 0; k < K; ++k)
            {
                const ContactJacEntry& entry = jacCache[k];

                OutDeriv vcA{};
                for (int b = 0; b < 2; ++b)
                {
                    const JacBlock& blk = entry.beam1Blocks[b];
                    vcA += (vel1[blk.frameIdx].getLinear()
                        + sofa::type::cross(vel1[blk.frameIdx].getAngular(),
                            blk.arm))
                        * blk.weight;
                }

                OutDeriv vcB{};
                for (int b = 0; b < 2; ++b)
                {
                    const JacBlock& blk = entry.beam2Blocks[b];
                    vcB += (vel2[blk.frameIdx].getLinear()
                        + sofa::type::cross(vel2[blk.frameIdx].getAngular(),
                            blk.arm))
                        * blk.weight;
                }
                
                const Vec3 dv = vcB - vcA;
                outVel[k] = Vec3(
                    s*(dv * entry.normal),
                    dv * entry.tangent1,
                    dv * entry.tangent2);
            }
        }
        else
        {
            // ── contactPoints mode: single interleaved output ─────────────────────
            //    outVel[0][2k]   = Ṗc_A[k]
            //    outVel[0][2k+1] = Ṗc_B[k]
            sofa::helper::WriteOnlyAccessor<sofa::core::objectmodel::Data<OutVecDeriv>>
                outVel = *dataVecOutVel[0];

            const sofa::Size newK2 = 2 * K;
            outVel.resize(newK2);
            for (sofa::Size k = 0; k < K; ++k)
            {
                const ContactJacEntry& entry = jacCache[k];

                OutDeriv vcA{};
                for (int b = 0; b < 2; ++b)
                {
                    const JacBlock& blk = entry.beam1Blocks[b];
                    vcA += (vel1[blk.frameIdx].getLinear()
                        + sofa::type::cross(vel1[blk.frameIdx].getAngular(),
                            blk.arm))
                        * blk.weight;
                }

                OutDeriv vcB{};
                for (int b = 0; b < 2; ++b)
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
    //  mode = "contactPoints":                                                  
    //    Single interleaved input force vector (size 2K):                      
    //    inForce[0][2k]   = FA at Pc_A[k]  →  Beam-1 frames:                   
    //      f += w_b·FA,   τ += w_b·(arm_A × FA)                                 
    //    inForce[0][2k+1] = FB at Pc_B[k]  →  Beam-2 frames:                  
    //      f += w_b·FB,   τ += w_b·(arm_B × FB)                                 
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

        if (!requireFrozenJacobianCache("applyJT(VecDeriv)"))
            return;

        sofa::helper::WriteAccessor<sofa::core::objectmodel::Data<In1VecDeriv>>
            outF1 = *dataVecOut1Force[0];
        sofa::helper::WriteAccessor<sofa::core::objectmodel::Data<In2VecDeriv>>
            outF2 = *dataVecOut2Force[0];

        const sofa::Size  K = static_cast<int>(m_jacCache.size());
        const bool gapMode = isGapMode();
        const Real s = d_gapSign.getValue();

        if (gapMode)
        {
            // Gap mode: single force input; back-project to both beams.
            sofa::helper::ReadAccessor<sofa::core::objectmodel::Data<OutVecDeriv>>
                inForce = *dataVecInForce[0];

            for (sofa::Size k = 0; k < K; ++k)
            {
                if (k >= inForce.size()) break;

                const ContactJacEntry& entry = m_jacCache[k];
                
                const Vec3 F = s * entry.normal * inForce[k][0]
                             + entry.tangent1 * inForce[k][1]
                             + entry.tangent2 * inForce[k][2];

                // Beam-1 receives negative contribution (Ṗc_A has − sign in δ̇).
                for (int b = 0; b < 2; ++b)
                {
                    const JacBlock& blk = entry.beam1Blocks[b];
                    outF1[blk.frameIdx].getLinear()
                        -= F * blk.weight;
                    outF1[blk.frameIdx].getAngular()
                        -= sofa::type::cross(blk.arm, F) * blk.weight;
                }

                // Beam-2 receives positive contribution.
                for (int b = 0; b < 2; ++b)
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

            for (sofa::Size k = 0; k < K; ++k)
            {
                const ContactJacEntry& entry = m_jacCache[k];

                // Force at Pc_A[k]: even slot → Beam-1.
                const sofa::Size slotA = 2 * k;
                if (slotA < inForce.size())
                {
                    const Vec3 FA = inForce[slotA];
                    for (int b = 0; b < 2; ++b)
                    {
                        const JacBlock& blk = entry.beam1Blocks[b];
                        outF1[blk.frameIdx].getLinear()
                            += FA * blk.weight;
                        outF1[blk.frameIdx].getAngular()
                            += sofa::type::cross(blk.arm, FA) * blk.weight;
                    }
                }

                // Force at Pc_B[k]: odd slot → Beam-2.
                const sofa::Size slotB = 2 * k + 1;
                if (slotB < inForce.size())
                {
                    const Vec3 FB = inForce[slotB];
                    for (int b = 0; b < 2; ++b)
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

    // ─────────────────────────────────────────────────────────────────────────
    //  applyJT  (MatrixDeriv)  –  constraint Jacobian assembly
    //
    //  Called by GenericConstraintSolver during the free-motion phase to build
    //  the per-constraint Jacobian rows passed to the LCP/QP solver.
    //
    //  mode = "contactPoints":
    //    Single input constraint matrix (from the one output MO).
    //    Column index encoding:  col = 2k   → Pc_A[k] constraint → Beam-1 only.
    //                            col = 2k+1 → Pc_B[k] constraint → Beam-2 only.
    //    In1 frame b (Beam-1): w_b·d (translational), w_b·(arm_A × d) (rotational).
    //    In2 frame b (Beam-2): w_b·d (translational), w_b·(arm_B × d) (rotational).
    //
    //  mode = "gap":                                                     
    //    dataMatIn[0] holds gap constraint rows.
    //    col k = contact index; d = constraint direction in contact-local frame.
    //    Back-project to world frame:
    //      d_phys = n̂·d[0] + t̂₁·d[1] + t̂₂·d[2]
    //    Then:
    //      In1 frame b: −w_b·d_phys (translational), −w_b·(arm_A × d_phys) (rotational).
    //      In2 frame b: +w_b·d_phys (translational), +w_b·(arm_B × d_phys) (rotational).
    //
    //    For DistanceUnilateralLagrangianConstraint, d = Vec3(1,0,0)
    //    → d_phys = n̂ (same as the old single-component implementation).
    //    Friction constraints writing Vec3(0,1,0) / Vec3(0,0,1) get t̂₁ / t̂₂.
    //
    //    Previously: d_phys = entry.normal * colIt.val()[0]  — only [0] used.
    // ─────────────────────────────────────────────────────────────────────────
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

        if (!requireFrozenJacobianCache("applyJT(MatrixDeriv)"))
            return;

        sofa::helper::WriteAccessor<sofa::core::objectmodel::Data<In1MatrixDeriv>>
            outM1 = *dataMatOut1[0];
        sofa::helper::WriteAccessor<sofa::core::objectmodel::Data<In2MatrixDeriv>>
            outM2 = *dataMatOut2[0];

        const Real s = d_gapSign.getValue();
    
        const sofa::Size  K = static_cast<int>(m_jacCache.size());
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
                    const sofa::Size k = colIt.index();
                    if (k < 0 || k >= K) continue;

                    const ContactJacEntry& entry = m_jacCache[k];
                    
                    const Vec3 d = s * entry.normal   * colIt.val()[0]
                                + entry.tangent1 * colIt.val()[1]
                                + entry.tangent2 * colIt.val()[2];

                    // Beam-1: negative Jacobian block (Ṗc_A enters with − in δ̇).
                    for (int b = 0; b < 2; ++b)
                    {
                        const JacBlock& blk = entry.beam1Blocks[b];
                        In1Deriv contrib;
                        contrib.getLinear() = -d * blk.weight;
                        contrib.getAngular() = -sofa::type::cross(blk.arm, d) * blk.weight;
                        row1.addCol(blk.frameIdx, contrib);
                    }

                    // Beam-2: positive Jacobian block.
                    for (int b = 0; b < 2; ++b)
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
                    const sofa::Size col = static_cast<sofa::Size>(colIt.index());
                    const sofa::Size k   = col / 2;
                    if (k < 0 || k >= K) continue;

                    const Vec3 d = colIt.val();
                    const ContactJacEntry& entry = m_jacCache[k];

                    if (col % 2 == 0)
                    {
                        // Even column → Pc_A[k] → Beam-1 frames only.
                        for (int b = 0; b < 2; ++b)
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
                        for (int b = 0; b < 2; ++b)
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
        factory->registerObjects(sofa::core::ObjectRegistrationData(
            "Maps SSIM contact-point descriptors (section IDs + curvilinear parameters) "
            "to two sets of Cosserat beam Rigid3d frames.\n"
            "Both modes use exactly ONE connected output MechanicalObject.\n"
            "Selectable via mappingMode:\n"
            "  'contactPoints': output size = 2K (interleaved). "
            "out[0][2k]=Pc_A (Beam-1 surface), out[0][2k+1]=Pc_B (Beam-2 surface). "
            "applyJT: even force/column indices -> Beam-1, odd -> Beam-2.\n"
            "  'gap': output size = K. "
            "applyJ gives gap velocity; applyJT back-projects with opposite signs to each beam.\n"
            "Implements apply / applyJ / applyJT(VecDeriv) / applyJT(MatrixDeriv) "
            "for use with FreeMotionAnimationLoop + GenericConstraintSolver.")
            .add<BeamContactMapping>());
    }

} // namespace Cosserat
