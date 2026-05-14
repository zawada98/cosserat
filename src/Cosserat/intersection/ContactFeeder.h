/******************************************************************************
 * Cosserat Plugin for SOFA Framework                                         *
 *                                                                            *
 * ContactFeeder.h                                                            *
 *                                                                            *
 * Bridges SphereSweptIntersectionMethod (SSIM) output to                    *
 * UnilateralLagrangianConstraint (or AugmentedLagrangianConstraint).        *
 *                                                                            *
 * At every AnimateBeginEvent this component:                                 *
 *   1. Reads the SSIM surface contact points and gap scalars.                *
 *   2. Filters pairs whose gap δn >= d_alarmDistance (no constraint needed). *
 *   3. Calls clear() + addContact() on the linked constraint object.         *
 *                                                                            *
 * Sign convention                                                            *
 * ---------------                                                            *
 *   norm[k]  = normalize( surfacePoints2[k] − surfacePoints1[k] )            *
 *   BaseContactLagrangianConstraint writes:                                  *
 *     c1[k][k] = −norm  (object1 = contactMO1, Beam-1 side)                *
 *     c2[k][k] = +norm  (object2 = contactMO2, Beam-2 side)                *
 *   getPositionViolation:                                                    *
 *     dfree = dot(Pc_B_free − Pc_A_free, norm) − contactDistance            *
 *     dfree < 0  ⟹  penetration  ⟹  λ ≥ 0 enforced                        *
 ******************************************************************************/
#pragma once

 // ── SOFA base infrastructure ──────────────────────────────────────────────────
 // BaseObject pulls in Data<T>, initData, msg_*, f_listening.
#include <sofa/core/objectmodel/BaseObject.h>
// Data.h is pulled in by BaseObject.h, but include explicitly to be safe.
#include <sofa/core/objectmodel/Data.h>
// SingleLink, initLink, BaseLink::FLAG_STRONGLINK.
#include <sofa/core/objectmodel/Link.h>
// Event base class (argument of handleEvent).
#include <sofa/core/objectmodel/Event.h>
// ObjectFactory forward declaration only; full header is in the .cpp.
namespace sofa::core { class ObjectFactory; }

// ── SOFA type library ─────────────────────────────────────────────────────────
#include <sofa/type/Vec.h>       // Vec3d
#include <sofa/type/vector.h>    // sofa::type::vector<T>

// ── SOFA default types ────────────────────────────────────────────────────────
#include <sofa/defaulttype/VecTypes.h>   // Vec3Types (needed for template param)

// ── Constraint type we populate ───────────────────────────────────────────────
#include <sofa/component/constraint/lagrangian/model/UnilateralLagrangianConstraint.h>

#include <Cosserat/engine/ContactTriad.h>

namespace Cosserat
{

    /**
     * @brief ContactFeeder
     *
     * Lightweight BaseObject that populates a
     * UnilateralLagrangianConstraint<Vec3Types> with per-contact data derived
     * from SphereSweptIntersectionMethod at the start of every time step.
     *
     * Usage (Python scene)
     * --------------------
     *   contact_node.addObject('ContactFeeder',
     *       surfacePoints1  = '@ssim.surfacePoints1',
     *       surfacePoints2  = '@ssim.surfacePoints2',
     *       distances       = '@ssim.distances',
     *       constraint      = '@contact_constraint',
     *       alarmDistance   = 0.005,
     *       mu              = 0.0)
     */
    class ContactFeeder : public sofa::core::objectmodel::BaseObject
    {
    public:
        SOFA_CLASS(ContactFeeder, sofa::core::objectmodel::BaseObject);

        // ── Convenience aliases ───────────────────────────────────────────────────
        using Vec3 = sofa::type::Vec3d;
        using Real = double;
        using VecVec3 = sofa::type::vector<Vec3>;

        using UnilateralConstraint =
            sofa::component::constraint::lagrangian::model
            ::UnilateralLagrangianConstraint<sofa::defaulttype::Vec3Types>;

        // Bring Data into unqualified scope inside this class so member
        // declarations stay readable. Fully qualified in the .cpp.
        template<class T>
        using Data = sofa::core::objectmodel::Data<T>;

        // ── Data inputs (linked from SSIM outputs) ────────────────────────────────
        
        /// Gap vectors {δn, 0, 0} per contact pair — only component [0] is used.
        /// δn = ‖Pint1 − Pint2‖ − (r1+r2).  δn < 0 ⇒ interpenetration.
        /// Link to SphereSweptIntersectionMethod::d_distances.
        Data<VecVec3> d_distances;

        // ── Parameters ────────────────────────────────────────────────────────────

        /// A contact constraint is generated for pair k when δn[k] < alarmDistance.
        /// Use a small positive value to catch near-contact before penetration.
        /// Default = 0.0  (only penetrating pairs generate constraints).
        Data<Real> d_alarmDistance;

        /// Coulomb friction coefficient μ.
        /// 0 = frictionless (1 DOF, UnilateralConstraintResolution).
        /// μ > 0 activates tangential friction (3 DOF, ...WithFriction).
        Data<Real> d_mu;
        
        /// Per-pair contact triad (n̂, t̂₁, t̂₂) produced by BCM.
        /// Read as  triads[k].n ,  triads[k].t1 ,  triads[k].t2 .
        /// Link to  @BCM.contactTriads .
        sofa::core::objectmodel::Data<VecContactTriad> d_contactTriads;

        /// Global gap sign s ∈ {+1, −1} such that (Pc_B − Pc_A)·n̂ = s·δn.
        /// Link to  @BCM.gapSign .
        sofa::core::objectmodel::Data<Real> d_gapSign;
        
        sofa::core::objectmodel::Data<Real> d_contactDistance;

        // ── Object link ───────────────────────────────────────────────────────────

        /// Link to the UnilateralLagrangianConstraint this feeder populates.
        /// object1 of that constraint must be contactMO1 (Beam-1 surface MO),
        /// object2 must be contactMO2 (Beam-2 surface MO).
        sofa::core::objectmodel::SingleLink<
            ContactFeeder,
            UnilateralConstraint,
            sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>  l_constraint;

        // ── SOFA lifecycle ────────────────────────────────────────────────────────
        ContactFeeder();
        ~ContactFeeder() override = default;

        void init()   override;
        void reinit() override;

        /// Intercepts AnimateBeginEvent to repopulate the constraint each step.
        void handleEvent(sofa::core::objectmodel::Event* e) override;

    private:
        /// Core routine – clear the constraint then add one entry per active pair.
        void feedContacts();

        static constexpr Real s_eps = Real(1e-14);
    };

    /// Factory registration – call from your plugin's initExternalModule().
    void registerContactFeeder(sofa::core::ObjectFactory* factory);

} // namespace Cosserat