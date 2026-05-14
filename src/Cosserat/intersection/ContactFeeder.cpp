/******************************************************************************
 * Cosserat Plugin for SOFA Framework                                         *
 *                                                                            *
 * ContactFeeder.cpp                                                          *
 *                                                                            *
 * See ContactFeeder.h for full documentation.                               *
 ******************************************************************************/
#include "ContactFeeder.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/simulation/AnimateBeginEvent.h>

#include <sofa/component/constraint/lagrangian/model/UnilateralConstraintResolution.h>

namespace Cosserat
{

// ─────────────────────────────────────────────────────────────────────────────
//  Constructor
// ─────────────────────────────────────────────────────────────────────────────
ContactFeeder::ContactFeeder()
    : Inherit1()
    , d_distances(initData(&d_distances,
        "distances",
        "Gap vectors {delta_n, 0, 0} per contact pair.\n"
        "Link to SphereSweptIntersectionMethod::d_distances.\n"
        "delta_n < 0 means interpenetration."))
    , d_alarmDistance(initData(&d_alarmDistance, Real(0.0),
        "alarmDistance",
        "Generate a contact constraint for pair k when delta_n[k] < alarmDistance.\n"
        "Use a small positive value (e.g. half an element length) to catch\n"
        "near-contacts before actual penetration occurs. Default = 0.0."))
    , d_mu(initData(&d_mu, Real(0.0),
        "mu",
        "Coulomb friction coefficient. 0 = frictionless (1 DOF per contact,\n"
        "UnilateralConstraintResolution). mu > 0 activates tangential friction\n"
        "(3 DOF per contact, UnilateralConstraintResolutionWithFriction)."))
    , l_constraint(initLink("constraint",
        "Link to the UnilateralLagrangianConstraint<Vec3Types> this feeder\n"
        "populates at each time step.\n"
        "Both object1 and object2 of that constraint must be the single\n"
        "interleaved BCM output MechanicalObject (size 2K).\n"
        "DOF 2k = Pc_A[k] (Beam-1 surface), DOF 2k+1 = Pc_B[k] (Beam-2 surface)."))
    , d_contactTriads(initData(&d_contactTriads,
        "contactTriads",
        "Per-pair contact triad (n̂, t̂₁, t̂₂) read as triads[k].n / .t1 / .t2.\n"
        "Link to @BCM.contactTriads."))
    , d_gapSign(initData(&d_gapSign, Real(1),
        "gapSign",
        "Global gap sign s ∈ {+1, −1} such that (Pc_B − Pc_A)·n̂ = s·δn.\n"
        "Link to @BCM.gapSign."))
    , d_contactDistance(initData(&d_contactDistance, Real(1e-3),
        "contactDistance",
        "Distance below which a contact is created"))
{
    // Enable event handling so handleEvent() is called by the simulation loop.
    this->f_listening.setValue(true);
}

// ─────────────────────────────────────────────────────────────────────────────
//  init
// ─────────────────────────────────────────────────────────────────────────────
void ContactFeeder::init()
{
    Inherit1::init();

    if (!l_constraint.get())
        msg_error() << "No UnilateralLagrangianConstraint linked via 'constraint'. "
                       "ContactFeeder will have no effect.";

    if (d_mu.getValue() < Real(0))
    {
        msg_warning() << "Friction coefficient mu = " << d_mu.getValue()
                      << " is negative. Clamping to 0 (frictionless).";
        d_mu.setValue(Real(0));
    }
	std::cout << "INITIALIZED SUCCESFFULLY";
}

// ─────────────────────────────────────────────────────────────────────────────
//  reinit
// ─────────────────────────────────────────────────────────────────────────────
void ContactFeeder::reinit()
{
    init();
}

// ─────────────────────────────────────────────────────────────────────────────
//  handleEvent
//
//  AnimateBeginEvent fires at the very start of each time step, before
//  FreeMotionAnimationLoop runs computeFreeMotion.  Repopulating the
//  constraint here guarantees that buildConstraintMatrix() and
//  getConstraintViolation() see the fresh contact list for this step.
// ─────────────────────────────────────────────────────────────────────────────
void ContactFeeder::handleEvent(sofa::core::objectmodel::Event* e)
{
    if (sofa::simulation::AnimateBeginEvent::checkEventType(e))
        feedContacts();
}

// ─────────────────────────────────────────────────────────────────────────────
//  feedContacts  –  core routine
//
//  For every contact pair k:
//
//    norm[k] = normalize( Pc_B[k] − Pc_A[k] )
//
//  This is the unit vector pointing from the Beam-1 surface toward the
//  Beam-2 surface.  BaseContactLagrangianConstraint stores it and uses it
//  in buildConstraintMatrix as:
//
//    c1[k][m1] += −norm       (object1 side, DOF m1 = 2k   in single MO)
//    c2[k][m2] += +norm       (object2 side, DOF m2 = 2k+1 in single MO)
//
//  and in getPositionViolation as:
//
//    dfree = dot( singleMO_free[2k+1] − singleMO_free[2k], norm ) − contactDistance
//          = dot( Pc_B_free − Pc_A_free, norm )
//
//  dfree < 0  ⇒  penetration  ⇒  UnilateralConstraintResolution clamps λ ≥ 0
//  dfree ≥ 0  ⇒  separation   ⇒  λ = 0 (no force)
//
//  BCM output MO layout (contactPoints mode, size 2K):
//    DOF 2k   = Pc_A[k]   (even index → Beam-1 surface)
//    DOF 2k+1 = Pc_B[k]   (odd  index → Beam-2 surface)
//
//  The MatrixDeriv column indices 2k and 2k+1 are then routed by
//  BeamContactMapping::applyJT(MatrixDeriv):
//    even col 2k   → Beam-1 Jacobian (+d, where d = −norm from ULC)
//    odd  col 2k+1 → Beam-2 Jacobian (+d, where d = +norm from ULC)
// ─────────────────────────────────────────────────────────────────────────────
void ContactFeeder::feedContacts()
{
    auto* cstr = l_constraint.get();
    if (!cstr) return;

    // Reading these Data triggers SSIM::doUpdate() / BCM::apply() if dirty.
    const auto& dists  = d_distances.getValue();

    // modified: read triad and gapSign from BCM instead of reconstructing the
    // normal from centerlines. The BCM-published normal is the geometric
    // ground truth for both external and nested contacts; gapSign is the
    // scalar that converts (Pc_B − Pc_A)·n̂_pub into the signed gap δn.
    // was: const auto& cl1 = d_centerlinePoints1.getValue();
    //      const auto& cl2 = d_centerlinePoints2.getValue();
    //      Vec3 n = cl2[k] - cl1[k];  ...  n = Vec3(0,0,1);
    const auto& triads = d_contactTriads.getValue();
    const Real  s_gap  = d_gapSign.getValue();
    

    const int K = static_cast<int>(triads.size());

    cstr->clear(K);
    if (K == 0) return;

    // modified: hard guard against an unlinked / wrongly-typed gapSign Data.
    // BCM must publish ±1 exactly; anything else silently corrupts the
    // constraint sign convention.
    if (std::abs(std::abs(s_gap) - Real(1)) > Real(1e-6))
    {
        msg_warning() << "gapSign = " << s_gap
                      << " is not ±1. ContactFeeder is most likely not linked "
                       "to @BCM.gapSign. Skipping all contacts this step.";
        return;
    }

    const Real alarm = d_alarmDistance.getValue();
    const Real mu    = d_mu.getValue();

    using Params = sofa::component::constraint::lagrangian::model
                   ::UnilateralLagrangianContactParameters;

    for (int k = 0; k < K; ++k)
    {
        const Real delta_n = dists[k][0];
        if (delta_n >= alarm) continue;     // separated more than alarm: no constraint

        // BCM convention:
        //   n̂_pub  = (Pc_B − Pc_A)/‖…‖           (always Beam-1 → Beam-2)
        //   (Pc_B − Pc_A) · n̂_pub  =  s · δn      (s = ±1, supplied by BCM)
        //
        // ULC computes  dfree = (xfree[m2] − xfree[m1]) · norm
        //                     = (Pc_B − Pc_A) · norm
        // Choosing norm = s · n̂_pub  ⇒  dfree = s² · δn = δn
        //   external (s=+1): norm = +n̂_pub
        //   nested   (s=−1): norm = −n̂_pub
        // In both cases dfree < 0 ⟺ penetration ⟹ λ ≥ 0 is enforced.

        const Vec3& n_pub = triads[k].n;
        const Real  n_len = n_pub.norm();
        if (n_len < s_eps) continue;       // degenerate triad — skip silently

        // modified: norm = s · n̂_pub. No local recomputation, no hardcoded axis.
        // was: n = cl2[k] - cl1[k]; n.normalize(); n = Vec3(0,0,1);
        const Vec3 norm = s_gap * n_pub;

        Params params(mu);
        cstr->addContact(params, norm,
                         /*contactDistance=*/Real(d_contactDistance.getValue()),
                         /*m1=*/2*k,
                         /*m2=*/2*k + 1);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  SOFA factory registration
// ─────────────────────────────────────────────────────────────────────────────
void registerContactFeeder(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(
        sofa::core::ObjectRegistrationData(
            "Populates a UnilateralLagrangianConstraint<Vec3Types> with per-contact\n"
            "data (normal direction, contact index, friction) from\n"
            "SphereSweptIntersectionMethod at the start of every time step.\n"
            // MODIFIED: description updated for single-MO interleaved layout.
            "Works with BeamContactMapping (contactPoints mode, single interleaved MO):\n"
            "  singleMO[2k] = Pc_A[k] (Beam-1 surface),\n"
            "  singleMO[2k+1] = Pc_B[k] (Beam-2 surface).\n"
            "The linked ULC must have object1 = object2 = singleMO.\n"
            "Full pipeline:\n"
            "  SSIM → ContactFeeder → UnilateralLagrangianConstraint\n"
            "       → BeamContactMapping::applyJT → Rigid3d frames\n"
            "       → GenericConstraintSolver")
        .add<ContactFeeder>());
}

} // namespace Cosserat