/******************************************************************************
 * Cosserat Plugin for SOFA Framework                                         *
 *                                                                            *
 * ContactPointsUnilateralConstraint.cpp                                      *
 *                                                                            *
 * See ContactPointsUnilateralConstraint.h for full documentation.           *
 ******************************************************************************/
#include "ContactPointsUnilateralConstraint.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/accessor.h>
#include <sofa/helper/rmath.h>
#include <sofa/type/RGBAColor.h>
#include <sofa/component/constraint/lagrangian/model/UnilateralConstraintResolution.h>

namespace Cosserat
{

using sofa::component::constraint::lagrangian::model::UnilateralConstraintResolution;
using sofa::component::constraint::lagrangian::model::UnilateralConstraintResolutionWithFriction;

// ─────────────────────────────────────────────────────────────────────────────
//  Constructor
// ─────────────────────────────────────────────────────────────────────────────
ContactPointsUnilateralConstraint::ContactPointsUnilateralConstraint()
    : Inherit()
    , d_contactTriads(initData(&d_contactTriads,
        "contactTriads",
        "Per-pair contact triad (n̂, t̂₁, t̂₂) read as triads[k].n / .t1 / .t2.\n"
        "Link to @BCM.contactTriads."))
    , d_gapSign(initData(&d_gapSign, Real(1),
        "gapSign",
        "Global gap sign s ∈ {+1, −1} such that (Pc_B − Pc_A)·n̂ = s·δn.\n"
        "Link to @BCM.gapSign."))
    , d_mu(initData(&d_mu, Real(0),
        "mu",
        "Coulomb friction coefficient μ.  0 → frictionless (1 row per pair);\n"
        "> 0 → friction (3 rows per pair, n̂ + t̂₁ + t̂₂)."))
{
}

// ─────────────────────────────────────────────────────────────────────────────
//  init / reinit
// ─────────────────────────────────────────────────────────────────────────────
void ContactPointsUnilateralConstraint::init()
{
    Inherit::init();

    if (!this->mstate)
    {
        msg_error() << "No MechanicalState<Vec3Types> resolved in the current "
                       "node.  Place this component in the same node as the "
                       "BeamContactMapping output contactMO.";
        return;
    }

    if (!d_contactTriads.getParent() && d_contactTriads.getValue().empty())
    {
        msg_warning() << "d_contactTriads is neither linked nor populated. "
                         "Link it to '@<bcm>.contactTriads' in the scene.";
    }

    if (!d_gapSign.getParent())
    {
        msg_warning() << "d_gapSign is not linked.  Link it to '@<bcm>.gapSign' "
                         "in the scene.  Using default " << d_gapSign.getValue()
                      << " for now.";
    }

    if (d_mu.getValue() < Real(0))
    {
        msg_warning() << "Friction coefficient mu = " << d_mu.getValue()
                      << " is negative.  Clamping to 0 (frictionless).";
        d_mu.setValue(Real(0));
    }
}

void ContactPointsUnilateralConstraint::reinit()
{
    init();
}

// ─────────────────────────────────────────────────────────────────────────────
//  rebuildContacts
//
//  Copies the per-pair triad (n, t1, t2) from d_contactTriads and caches it
//  per contact.  Zero arithmetic between producer (BCM) and consumer (CPULC).
void ContactPointsUnilateralConstraint::rebuildContacts()
{
    m_constraintRowsBuilt = false;
    const auto& triads = d_contactTriads.getValue();
    const size_t K     = triads.size();

    m_contacts.clear();
    m_sign = d_gapSign.getValue();

    if (K == 0) return;

    // Validate gap sign is within a small tolerance of ±1.
    if (std::abs(std::abs(m_sign) - Real(1)) > Real(1e-6))
    {
        msg_warning() << "d_gapSign = " << m_sign
                      << " is not ±1.  Snapping to sign.";
        m_sign = (m_sign >= Real(0)) ? Real(+1) : Real(-1);
    }

    // was `if (moSize != 2 * K) { msg_error(); return; }`.
    // Strict equality is brittle: it fires whenever BCM's triad vector and
    // its output MO get out of sync by any amount, even when every pair k we
    // care about still has valid slots 2k and 2k+1.  The real per-pair
    // invariant is checked inside the loop below.
    const size_t moSize = this->mstate->getSize();                       
    const bool withFriction = d_mu.getValue() > Real(0);
    m_contacts.reserve(K);

    for (size_t k = 0; k < K; ++k)
    {
        //per-pair index guard.  Each pair needs the two
        // interleaved slots 2k (Pc_A) and 2k+1 (Pc_B) to exist in contactMO.
        // BCM resizes contactMO monotonically, so a miss here means BCM did
        // not run this step, or CPULC is wired to the wrong MO.  Subsequent
        // pairs have strictly larger indices and would all fail the same
        // way, so `break` is sufficient and avoids error spam.
        if (2 * k + 1 >= moSize)                                        
        {                                                                 
            msg_error() << "Contact pair " << k                           
                        << " requires contactMO slots " << (2 * k)        
                        << " and " << (2 * k + 1)                       
                        << " but mstate->getSize()=" << moSize << ".  "  
                        << "Either BCM did not run this step, or CPULC " 
                        << "is linked to a different MO than BCM's "     
                        << "contactPoints output.  Registering only "   
                        << "pairs [0," << k << ").";                      
            break;                                                     
        }                                                                

        const ContactTriad& triad = triads[k];

        // BCM publishes a zero triad for degenerate pairs.  Skip them.
        if (triad.n.norm2() < Real(1e-24))
            continue;

        Contact c;
        c.k       = static_cast<int>(k);
        c.n       = triad.n;
        c.cId     = 0;
        c.dfree_n = Real(0);

        if (withFriction)
        {
            c.t1 = triad.t1;
            c.t2 = triad.t2;
        }

        m_contacts.push_back(c);
    }
}
// ─────────────────────────────────────────────────────────────────────────────
//  buildConstraintMatrix
// ─────────────────────────────────────────────────────────────────────────────
void ContactPointsUnilateralConstraint::buildConstraintMatrix(
    const sofa::core::ConstraintParams* /*cParams*/,
    DataMatrixDeriv& c_d,
    unsigned int&    cIndex,
    const DataVecCoord& /*x*/)
{
    rebuildContacts();
    
    if (m_contacts.empty()) return;

    const bool withFriction = d_mu.getValue() > Real(0);
    auto c = sofa::helper::getWriteAccessor(c_d);

    for (Contact& con : m_contacts)
    {
        const int  m1 = 2 * con.k;
        const int  m2 = 2 * con.k + 1;
        const Vec3 sn = m_sign * con.n;

        con.cId = cIndex++;

        auto row_n = c->writeLine(con.cId);
        row_n.addCol(m1, -sn);
        row_n.addCol(m2,  sn);

        if (withFriction)
        {
            auto row_t = c->writeLine(con.cId + 1);
            row_t.addCol(m1, -con.t1);
            row_t.addCol(m2,  con.t1);

            auto row_s = c->writeLine(con.cId + 2);
            row_s.addCol(m1, -con.t2);
            row_s.addCol(m2,  con.t2);

            cIndex += 2;
        }
    }
    m_constraintRowsBuilt = true;   
}

// ─────────────────────────────────────────────────────────────────────────────
//  getConstraintViolation
// ─────────────────────────────────────────────────────────────────────────────

void ContactPointsUnilateralConstraint::getConstraintViolation(
    const sofa::core::ConstraintParams* cParams,
    sofa::linearalgebra::BaseVector*    v,
    const DataVecCoord&                 xfree_d,
    const DataVecDeriv&                 vfree_d)
{
    if (!cParams) return;                 
    if (!m_constraintRowsBuilt) return;    
    if (m_contacts.empty()) return;

    const VecCoord& xfree = xfree_d.getValue();
    const VecCoord& x     = this->mstate->read(
                                sofa::core::vec_id::read_access::position)->getValue();

    const bool posOrder =
        (cParams->constOrder() == sofa::core::ConstraintOrder::POS ||
         cParams->constOrder() == sofa::core::ConstraintOrder::POS_AND_VEL);

    const bool withFriction = d_mu.getValue() > Real(0);

    if (posOrder)
    {
        for (Contact& con : m_contacts)
        {
            const int m1 = 2 * con.k;
            const int m2 = 2 * con.k + 1;

            const Coord& Pfree = xfree[m2];
            const Coord& Qfree = xfree[m1];
            const Coord& P     = x    [m2];
            const Coord& Q     = x    [m1];

            const Coord PPfree = Pfree - P;
            const Coord QQfree = Qfree - Q;
            const Real  ref_dist = PPfree.norm() + QQfree.norm();

            const Real dfree = m_sign * sofa::type::dot(Pfree - Qfree, con.n);
            const Real delta = m_sign * sofa::type::dot(P     - Q    , con.n);

            v->set(con.cId, dfree);
            con.dfree_n = dfree;

            if (!withFriction) continue;

            // All tangential violations below use bare con.t1/ con.t2 —    
            // no m_sign factor, matching the tangent rows in             
            // buildConstraintMatrix.                                        
            Real dfree_t = Real(0), dfree_s = Real(0);

            const bool bothNearZero =
                sofa::helper::rabs(delta) < Real(1e-5) * ref_dist &&
                sofa::helper::rabs(dfree) < Real(1e-5) * ref_dist;

            if (bothNearZero)
            {
                dfree_t = sofa::type::dot(PPfree, con.t1) -      
                          sofa::type::dot(QQfree, con.t1);
                dfree_s = sofa::type::dot(PPfree, con.t2) -       
                          sofa::type::dot(QQfree, con.t2);
            }
            else if (sofa::helper::rabs(delta - dfree) >
                     Real(1e-3) * sofa::helper::rabs(delta))
            {
                const Real dt_frac = delta / (delta - dfree);
                if (dt_frac > Real(0) && dt_frac < Real(1))
                {
                    const Coord Pt = P * (Real(1) - dt_frac) + Pfree * dt_frac;
                    const Coord Qt = Q * (Real(1) - dt_frac) + Qfree * dt_frac;
                    const Coord PtPfree = Pfree - Pt;
                    const Coord QtQfree = Qfree - Qt;
                    dfree_t = sofa::type::dot(PtPfree, con.t1) -  
                              sofa::type::dot(QtQfree, con.t1);
                    dfree_s = sofa::type::dot(PtPfree, con.t2) -  
                              sofa::type::dot(QtQfree, con.t2);
                }
                else if (dfree < Real(0))
                {
                    dfree_t = sofa::type::dot(PPfree, con.t1) - 
                              sofa::type::dot(QQfree, con.t1);
                    dfree_s = sofa::type::dot(PPfree, con.t2) -  
                              sofa::type::dot(QQfree, con.t2);
                }
                // else: contact released mid-step → tangentials stay 0.
            }
            else
            {
                if (dfree < Real(0))
                {
                    dfree_t = sofa::type::dot(PPfree, con.t1) -       
                              sofa::type::dot(QQfree, con.t1);
                    dfree_s = sofa::type::dot(PPfree, con.t2) -       
                              sofa::type::dot(QQfree, con.t2);
                }
            }

            v->set(con.cId + 1, dfree_t);
            v->set(con.cId + 2, dfree_s);
        }
    }
    else  // VEL or ACC order
    {
        const VecDeriv& vfree = vfree_d.getValue();
        const SReal     dt    = this->getContext()->getDt();
        const SReal     invDt = SReal(1) / dt;

        for (const Contact& con : m_contacts)
        {
            const int m1 = 2 * con.k;
            const int m2 = 2 * con.k + 1;

            const Deriv QP_invDt = (x[m2] - x[m1]) * invDt;
            const Deriv QP_vfree = vfree[m2] - vfree[m1];
            const Deriv dVec     = QP_vfree + QP_invDt;

            v->set(con.cId, m_sign * sofa::type::dot(dVec, con.n));

            if (withFriction)
            {
                v->set(con.cId + 1, sofa::type::dot(QP_vfree, con.t1)); 
                v->set(con.cId + 2, sofa::type::dot(QP_vfree, con.t2)); 
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  getConstraintResolution
// ─────────────────────────────────────────────────────────────────────────────
void ContactPointsUnilateralConstraint::getConstraintResolution(
    const sofa::core::ConstraintParams* /*cParams*/,
    std::vector<sofa::core::behavior::ConstraintResolution*>& resTab,
    unsigned int& offset)
{
    const Real mu           = d_mu.getValue();
    const bool withFriction = mu > Real(0);
    const unsigned int rowsPerPair  = withFriction ? 3u : 1u; 
    const size_t neededSize    = static_cast<size_t>(offset) +
                                 m_contacts.size() * rowsPerPair; 
    
    if (neededSize > resTab.size())                                
    {
        msg_error() << "getConstraintResolution: resTab size "
                    << resTab.size() << " < needed " << neededSize
                    << " (offset=" << offset << ", pairs="
                    << m_contacts.size() << ", rowsPerPair="
                    << rowsPerPair << "). Upstream constraint solver "
                    << "did not pre-size resTab. Aborting this "
                    << "component's registration to prevent heap "
                    << "corruption. Constraint will be silently "
                    << "ignored this step.";
        // Keep offset consistent with the row count written in
        // buildConstraintMatrix so the solver's bookkeeping doesn't drift.
        offset += static_cast<unsigned int>(m_contacts.size() * rowsPerPair);
        return;
    }

    for (const Contact& con : m_contacts)
    {
        if (withFriction)
        {
            resTab[offset] = new UnilateralConstraintResolutionWithFriction(mu);
            offset += 3;
        }
        else
        {
            resTab[offset] = new UnilateralConstraintResolution();
            offset += 1;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  isActive
// ─────────────────────────────────────────────────────────────────────────────
bool ContactPointsUnilateralConstraint::isActive() const
{
    return !m_contacts.empty();  
}

// ─────────────────────────────────────────────────────────────────────────────
//  draw — normal + (if friction) tangents, rooted at each surface point.
// ─────────────────────────────────────────────────────────────────────────────
void ContactPointsUnilateralConstraint::draw(
    const sofa::core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowInteractionForceFields()) return;
    if (m_contacts.empty() || !this->mstate) return;

    const VecCoord& x = this->mstate->read(
                            sofa::core::vec_id::read_access::position)->getValue();

    const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();
    vparams->drawTool()->disableLighting();

    const bool withFriction = d_mu.getValue() > Real(0);

    std::vector<sofa::type::Vec3>      verts;
    std::vector<sofa::type::RGBAColor> cols;

    for (const Contact& con : m_contacts)
    {
        const int m1 = 2 * con.k;
        const int m2 = 2 * con.k + 1;
        if (m1 >= static_cast<int>(x.size()) ||
            m2 >= static_cast<int>(x.size())) continue;

        const Vec3 Pa = x[m1];
        const Vec3 Pb = x[m2];

        // Normal at each surface, drawn outward then returning.
        verts.push_back(Pa); verts.push_back(Pa + con.n);
        cols.push_back(sofa::type::RGBAColor::white());
        verts.push_back(Pb); verts.push_back(Pb - con.n);
        cols.push_back(sofa::type::RGBAColor(0, 0.5, 0.5, 1));

        if (withFriction)
        {
            // Tangents drawn at the midpoint at half length.
            const Vec3 M = (Pa + Pb) * Real(0.5);
            verts.push_back(M); verts.push_back(M + con.t1 * Real(0.5));
            cols.push_back(sofa::type::RGBAColor::red());
            verts.push_back(M); verts.push_back(M + con.t2 * Real(0.5));
            cols.push_back(sofa::type::RGBAColor::green());
        }
    }

    vparams->drawTool()->drawLines(verts, 3, cols);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Factory registration
// ─────────────────────────────────────────────────────────────────────────────
void registerContactPointsUnilateralConstraint(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(
        sofa::core::ObjectRegistrationData(
            "Single-MO unilateral contact constraint operating on the\n"
            "BeamContactMapping contactPoints-mode output MO.\n"
            "\n"
            "DOF layout (interleaved contactMO, size 2K):\n"
            "  [2k]   = Pc_A[k]  (Beam-1 surface)\n"
            "  [2k+1] = Pc_B[k]  (Beam-2 surface)\n"
            "\n"
            "Inputs (linked from BCM):\n"
            "  contactTriads : vector<ContactTriad> — one (n, t1, t2) per pair.\n"
            "  gapSign       : Real, ±1 (makes dfree<0 ⟺ penetration).\n"
            "\n"
            "Constraint (per pair k):\n"
            "  dfree[k] = s · n̂[k] · (xfree[2k+1] − xfree[2k])\n"
            "  μ = 0  → 1 row per pair, UnilateralConstraintResolution.\n"
            "  μ > 0 → 3 rows per pair (n̂ + t̂₁ + t̂₂),\n"
            "          UnilateralConstraintResolutionWithFriction.\n"
            "\n"
            "Scene usage:\n"
            "  contact_node.addObject(\n"
            "      'ContactPointsUnilateralConstraint',\n"
            "      name          = 'cpulc',\n"
            "      contactTriads = '@../bcm.contactTriads',\n"
            "      gapSign       = '@../bcm.gapSign',\n"
            "      mu            = 0.0)")
        .add<ContactPointsUnilateralConstraint>());
}

} // namespace Cosserat
