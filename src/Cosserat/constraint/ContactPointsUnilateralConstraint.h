/******************************************************************************
 * Cosserat Plugin for SOFA Framework                                         *
 *                                                                            *
 * ContactPointsUnilateralConstraint.h                                        *
 *                                                                            *
 * Single-MO unilateral contact constraint for the BeamContactMapping         *
 * "contactPoints" pipeline.                                                  *
 *                                                                            *
 * Inherits from sofa::core::behavior::LagrangianConstraint<Vec3Types>        *
 * (single-MO), operating on the interleaved contactMO produced by BCM:       *
 *     out[0][2k]   = Pc_A[k]   (Beam-1 surface)                              *
 *     out[0][2k+1] = Pc_B[k]   (Beam-2 surface)                              *
 *                                                                            *
 * ── Constraint equation (per pair k) ─────────────────────────────────────── *
 *   dfree[k] = s · n̂[k] · (xfree[2k+1] − xfree[2k])                         *
 *   dfree[k] < 0  ⟺  penetration  ⟹  λ ≥ 0 via UnilateralConstraintResolution *
 *   μ > 0 adds two tangential rows per pair with Coulomb friction.           *
 *                                                                            *
 * ── Inputs from the scene (all linked from BCM) ─────────────────────────── *
 *                                                                            *
 *   contactTriads → @BCM.contactTriads               *
 *                                                                            *
 *     CPULC NEVER constructs tangents locally — the triad comes from BCM,    *
 *     which gets it from the SSIM contact geometry (well-defined directions).*
 *                                                                            *
 *   gapSign       → @BCM.gapSign         Real ∈ {+1, −1}                     *
 *     Global sign s such that (Pc_B − Pc_A)·n̂ = s·δn.  Makes dfree < 0 ⟺    *
 *     penetration regardless of which beam is external/internal.             *
 *                                                                            *
 *   mu            → Coulomb friction coefficient (Real, default 0).          *
 ******************************************************************************/
#pragma once
#include <Cosserat/config.h>
#include <Cosserat/engine/ContactTriad.h>
#include <sofa/core/behavior/LagrangianConstraint.h>
#include <sofa/core/behavior/ConstraintResolution.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>

namespace sofa::core { class ObjectFactory; }

namespace Cosserat
{

class SOFA_COSSERAT_API ContactPointsUnilateralConstraint
    : public sofa::core::behavior::LagrangianConstraint<sofa::defaulttype::Vec3Types>
{
public:
    SOFA_CLASS(
        ContactPointsUnilateralConstraint,
        SOFA_TEMPLATE(sofa::core::behavior::LagrangianConstraint,
                      sofa::defaulttype::Vec3Types));

    // ── Type aliases ─────────────────────────────────────────────────────────
    using Inherit         = sofa::core::behavior::LagrangianConstraint<
                                sofa::defaulttype::Vec3Types>;
    using DataTypes       = sofa::defaulttype::Vec3Types;
    using Real            = DataTypes::Real;
    using Coord           = DataTypes::Coord;
    using Deriv           = DataTypes::Deriv;
    using VecCoord        = DataTypes::VecCoord;
    using VecDeriv        = DataTypes::VecDeriv;
    using MatrixDeriv     = DataTypes::MatrixDeriv;
    using DataVecCoord    = Inherit::DataVecCoord;
    using DataVecDeriv    = Inherit::DataVecDeriv;
    using DataMatrixDeriv = Inherit::DataMatrixDeriv;
    using Vec3            = sofa::type::Vec3d;

    // ── Data inputs (linked from BCM) ────────────────────────────────────────

    /// Per-pair contact triad (n̂, t̂₁, t̂₂) produced by BCM.
    /// Read as  triads[k].n ,  triads[k].t1 ,  triads[k].t2 .
    /// Link to  @BCM.contactTriads .
    sofa::core::objectmodel::Data<VecContactTriad> d_contactTriads;

    /// Global gap sign s ∈ {+1, −1} such that (Pc_B − Pc_A)·n̂ = s·δn.
    /// Link to  @BCM.gapSign .
    sofa::core::objectmodel::Data<Real> d_gapSign;

    // ── Parameter ────────────────────────────────────────────────────────────

    /// Coulomb friction coefficient μ.
    ///   0  → frictionless → 1 constraint row per pair.
    ///   >0 → friction     → 3 constraint rows per pair.
    sofa::core::objectmodel::Data<Real> d_mu;

    // ── SOFA lifecycle ───────────────────────────────────────────────────────
    ContactPointsUnilateralConstraint();
    ~ContactPointsUnilateralConstraint() override = default;

    void init()   override;
    void reinit() override;

    // ── LagrangianConstraint overrides ───────────────────────────────────────

    void buildConstraintMatrix(
        const sofa::core::ConstraintParams* cParams,
        DataMatrixDeriv& c,
        unsigned int&    cIndex,
        const DataVecCoord& x) override;

    void getConstraintViolation(
        const sofa::core::ConstraintParams* cParams,
        sofa::linearalgebra::BaseVector* v,
        const DataVecCoord& x,
        const DataVecDeriv& v_vel) override;

    void getConstraintResolution(
        const sofa::core::ConstraintParams* cParams,
        std::vector<sofa::core::behavior::ConstraintResolution*>& resTab,
        unsigned int& offset) override;

    void draw(const sofa::core::visual::VisualParams* vparams) override;

    bool isActive() const override;

    sofa::type::vector<std::string> getConstraintIdentifiers() override
    {
        return { "Unilateral", "ContactPoints" };
    }

private:
    struct Contact
    {
        int          k;       ///< pair index in SSIM/BCM outputs
        Vec3         n;       ///< n̂[k]   
        Vec3         t1;       ///< t̂₁[k]   (valid only when μ > 0)
        Vec3         t2;       ///< t̂₂[k]   (valid only when μ > 0)
        unsigned int cId;     ///< constraint row id written by buildConstraintMatrix
        Real         dfree_n; ///< last computed normal gap, for isActive()
    };

    sofa::type::vector<Contact> m_contacts;

    /// Scalar gap sign cached once per buildConstraintMatrix call.
    Real m_sign = Real(1);

    /// Rebuilds m_contacts from d_contactTriads and d_gapSign.  Validates the
    /// MO size (mstate->getSize() must equal 2·K) and reports a hard error on
    /// mismatch.  Never constructs tangents locally — triad is extracted via
    /// Quat::rotate on the published frame orientation.
    void rebuildContacts();
    
    /// True once buildConstraintMatrix has run and assigned valid cIds.  
    /// Reset to false at the start of rebuildContacts (and therefore at    
    /// the top of every buildConstraintMatrix call). getConstraintViolation 
    /// bails out if this is false to avoid writing violations at stale cIds 
    /// held over from the previous step.                                    
    bool m_constraintRowsBuilt = false;                                       
};

/// Factory registration — call from the Cosserat plugin initExternalModule().
void registerContactPointsUnilateralConstraint(sofa::core::ObjectFactory* factory);

} // namespace Cosserat