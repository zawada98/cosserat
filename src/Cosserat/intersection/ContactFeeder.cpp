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
    , d_surfacePoints1(initData(&d_surfacePoints1,
        "surfacePoints1",
        "Surface contact points on Beam 1 (Pc_A[k]).\n"
        "Link to SphereSweptIntersectionMethod::d_surfacePoints1."))
    , d_surfacePoints2(initData(&d_surfacePoints2,
        "surfacePoints2",
        "Surface contact points on Beam 2 (Pc_B[k]).\n"
        "Link to SphereSweptIntersectionMethod::d_surfacePoints2."))
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
    , d_centerlinePoints1(initData(&d_centerlinePoints1,
    "centerlinePoints1",
    "Centreline contact points on Beam 1 (Pint_A[k]).\n"
    "Link to SphereSweptIntersectionMethod::d_centerlinePoints1.\n"
    "Used to compute the sign-safe contact normal."))
    , d_centerlinePoints2(initData(&d_centerlinePoints2,
    "centerlinePoints2",
    "Centreline contact points on Beam 2 (Pint_B[k]).\n"
    "Link to SphereSweptIntersectionMethod::d_centerlinePoints2.\n"
    "Used to compute the sign-safe contact normal."))
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

    // Reading these Data triggers SSIM::doUpdate() if dirty.
    const auto& pts1  = d_surfacePoints1.getValue();
    const auto& pts2  = d_surfacePoints2.getValue();
    const auto& cl1 = d_centerlinePoints1.getValue();
    const auto& cl2 = d_centerlinePoints2.getValue();
    const auto& dists = d_distances.getValue();

    const int K = static_cast<int>(
    std::min({ pts1.size(), pts2.size(), dists.size(),
               cl1.size(), cl2.size() }));

    // Clear the previous step's contacts and reserve space for up to K entries.
    cstr->clear(K);

    if (K == 0) return;

    const Real alarm = d_alarmDistance.getValue();
    const Real mu    = d_mu.getValue();

    using Params =
        sofa::component::constraint::lagrangian::model
            ::UnilateralLagrangianContactParameters;

    for (int k = 0; k < K; ++k)
    {
        const Real delta_n = dists[k][0];
        if (delta_n >= alarm) continue; //todo

        // Contact normal: unit vector from Pc_A (Beam-1 surface) to Pc_B (Beam-2 surface).
        // Falls back to skipping this pair if the surface points are coincident
        // (degenerate case; BeamContactMapping will also have logged a warning).
        Vec3 n = cl2[k] - cl1[k];
        const Real d = n.norm();
        if (d < s_eps) continue;
        n = Vec3(0, 0, 1);//todo

        // addContact(params, norm, contactDistance, m1, m2)
        //
        // MODIFIED: BCM now uses a single interleaved output MO (size 2K).
        //   m1 = 2k   → singleMO[2k]   = Pc_A[k]  (Beam-1 surface, object1 side)
        //   m2 = 2k+1 → singleMO[2k+1] = Pc_B[k]  (Beam-2 surface, object2 side)
        //
        // Previously: m1 = k, m2 = k   (with two separate MOs of size K each).
        Params params(mu);
        cstr->addContact(params, n,
            /*contactDistance=*/Real(0.0),
                         2*k,
                         2*k+1);
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