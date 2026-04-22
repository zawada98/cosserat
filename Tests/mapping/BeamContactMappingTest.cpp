/******************************************************************************
 * Cosserat Plugin for SOFA Framework
 *
 * BeamContactMappingTest.cpp
 *
 * C++ GTest suite for BeamContactMapping (BCM).
 *
 * ── No scene graph ───────────────────────────────────────────────────────────
 *
 * The previous version of this test built a full SOFA scene (nodes,
 * MechanicalObjects, simpleapi) just to put the BCM into a callable state.
 * That scaffolding was unnecessary.
 *
 * apply(), applyJ() and applyJT() are all PUBLIC methods whose signatures
 * explicitly receive their input and output Data objects as arguments:
 *
 *   apply   (mparams, dataVecOutPos,  dataVecIn1Pos, dataVecIn2Pos)
 *   applyJ  (mparams, dataVecOutVel,  dataVecIn1Vel, dataVecIn2Vel)
 *   applyJT (mparams, dataVecOut1Frc, dataVecOut2Frc, dataVecInFrc)
 *
 * We therefore create standalone Data<T> objects as local variables, fill
 * them with the test geometry, call apply() to populate m_jacCache, and
 * then call applyJ() / applyJT() directly.  No SOFA node, no simulation,
 * no scene graph is needed at any point.
 *
 * ── Why apply() must be called before applyJ/applyJT ────────────────────────
 *
 * applyJ() and applyJT() do not recompute geometry.  They read m_jacCache,
 * a private member that stores the per-contact Jacobian blocks (frame
 * indices, interpolation weights, moment arms) computed by apply().
 * Calling applyJ() with an empty m_jacCache produces zero output.
 * apply() must therefore always be called first.
 *
 * ── How expected values are obtained (no reimplementation) ──────────────────
 *
 * apply() — gap mode
 *   Tests use geometries where the expected gap is obvious by inspection:
 *
 *   Configuration A  (z-separation)
 *     Beam-1 frames at z=0, Beam-2 frames at z=H, same x,y on both beams.
 *     Any interpolated P_A and P_B are directly above each other → n̂=(0,0,1).
 *     By the DEFINITION of gap = surface-to-surface distance along n̂:
 *         gap = (0, 0, H − R1 − R2)
 *     No formula needed; this follows from what "gap" means.
 *
 *   Configuration B  (x-separation)
 *     Same argument, separation D along x → gap = (D − R1 − R2, 0, 0).
 *
 * apply() — contactPoints mode
 *   The same geometries are used, but the individual surface points are
 *   checked rather than their difference:
 *
 *       Pc_A = P_A + R1·n̂        (Beam-1 surface displaced toward Beam-2)
 *       Pc_B = P_B − R2·n̂        (Beam-2 surface displaced toward Beam-1)
 *
 *   For configuration A (n̂=(0,0,1)):
 *       Pc_A = (P_A.x, P_A.y, R1)
 *       Pc_B = (P_B.x, P_B.y, H − R2)
 *
 *   A cross-check verifies Pc_B − Pc_A == gap-mode output.
 *
 * applyJ() / applyJT()
 *   (i)  Symmetry:        max|J − JT.T| ≤ TOL_J
 *   (ii) Virtual work:    |F·(J·v) − v·(JT·F)| ≤ TOL_VW
 *   (iii)Translational spot-check — gap mode:
 *        δ̇ = Ṗc_B − Ṗc_A, so:
 *          J[gap_z, vz_beam1_i] = −(1−α),   J[gap_z, vz_beam2_j] = +(1−β).
 *   (iv) Translational spot-check — contactPoints mode (n̂ frozen):
 *        J layout: rows 0..3K-1 = Ṗc_A,  rows 3K..6K-1 = Ṗc_B.
 *        Pc_A depends only on Beam-1 → J_A[Pc_A.z, vz_beam1_i] = +(1−α).
 *        Pc_B depends only on Beam-2 → J_B[Pc_B.z, vz_beam2_j] = +(1−β).
 *        Cross-coupling must be zero in both rows.
 *   (v)  Gap-equals-difference: J_gap == J_cp_B − J_cp_A entry-by-entry.
 *   (vi) ALGO_2 sparsity: exactly one non-zero Beam-1 block per contact.
 *   (vii)contactPoints no-cross-beam coupling: Pc_A rows silent on Beam-2
 *        columns, and vice versa (frozen n̂ linearisation).
 ******************************************************************************/

#include <gtest/gtest.h>

#include <Cosserat/src/Cosserat/mapping/BeamContactMapping.h>

// SOFA object infrastructure — needed to instantiate BeamContactMapping
// without a scene graph.
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Quat.h>

#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <array>
#include <algorithm>

// ════════════════════════════════════════════════════════════════════════════
//  Sec 1 — Type aliases
// ════════════════════════════════════════════════════════════════════════════

namespace
{

using Rigid3dTypes = sofa::defaulttype::Rigid3dTypes;
using Vec3dTypes   = sofa::defaulttype::Vec3dTypes;

using In1VecCoord = Rigid3dTypes::VecCoord;
using In1VecDeriv = Rigid3dTypes::VecDeriv;  // In2Types == In1Types == Rigid3d
using In2VecDeriv = Rigid3dTypes::VecDeriv;
using In1Coord    = Rigid3dTypes::Coord;

using OutVecCoord = Vec3dTypes::VecCoord;
using OutVecDeriv = Vec3dTypes::VecDeriv;

using Vec3d = sofa::type::Vec3d;
using Vec2d = sofa::type::Vec2d;
using Vec2i = sofa::type::Vec<2, int>;

// Shorthand for the Data types the BCM methods expect.
using D_In1Coord = sofa::core::objectmodel::Data<In1VecCoord>;
using D_In1Deriv = sofa::core::objectmodel::Data<In1VecDeriv>;
using D_In2Deriv = sofa::core::objectmodel::Data<In2VecDeriv>;
using D_OutCoord = sofa::core::objectmodel::Data<OutVecCoord>;
using D_OutDeriv = sofa::core::objectmodel::Data<OutVecDeriv>;

// ════════════════════════════════════════════════════════════════════════════
//  Sec 2 — Test constants
// ════════════════════════════════════════════════════════════════════════════

constexpr double R1 = 0.10;  // Beam-1 cross-section radius
constexpr double R2 = 0.05;  // Beam-2 cross-section radius (≠ R1: catches sign bugs)
constexpr double H  = 0.50;  // z-separation  (H  > R1+R2 = 0.15 → separated)
constexpr double D  = 0.30;  // x-separation  (D  > R1+R2 = 0.15 → separated)

constexpr double TOL_APPLY = 1e-9;
constexpr double TOL_J     = 1e-9;
constexpr double TOL_VW    = 1e-9;

constexpr int    RNG_SEED  = 42;

// ════════════════════════════════════════════════════════════════════════════
//  Sec 3 — Geometry helpers  (scene setup only; no mapping logic)
// ════════════════════════════════════════════════════════════════════════════

In1Coord makeRigid(double x, double y, double z,
                   double qx=0.0, double qy=0.0, double qz=0.0, double qw=1.0)
{
    In1Coord c;
    c.getCenter()      = Vec3d(x, y, z);
    c.getOrientation() = sofa::type::Quat<double>(qx, qy, qz, qw);
    return c;
}

/// Quaternion for rotation by deg degrees about z.
std::array<double,4> rotZ(double deg)
{
    const double h = deg * M_PI / 360.0;
    return {0.0, 0.0, std::sin(h), std::cos(h)};
}

// ════════════════════════════════════════════════════════════════════════════
//  Sec 4 — BcmFixture
//
//  Owns a single BeamContactMapping instance and the standalone Data
//  objects that the mapping methods read from / write to.
//
//  Why no scene graph?
//  ───────────────────
//  apply(), applyJ() and applyJT() are public methods that take ALL their
//  inputs and outputs as explicit Data<T>* pointer-vector arguments.  They
//  do NOT internally query the scene graph during execution.
//
//  Mapping modes
//  ─────────────
//  "gap" (default):
//    apply()   writes one output MO:   d_out[k]  = δ[k] = Pc_B − Pc_A
//    applyJ()  reads one output slot:  d_vout
//    applyJT() reads one input slot:   d_fout
//
//    J shape:  (3K) × (6N1+6N2)
//    JT shape: (6N1+6N2) × (3K)
//
//  "contactPoints":
//    apply()   writes two output MOs:  d_out[k]  = Pc_A[k]
//                                      d_outB[k] = Pc_B[k]
//    applyJ()  reads two output slots: d_vout  (Ṗc_A, rows 0..3K-1)
//                                      d_voutB (Ṗc_B, rows 3K..6K-1)
//    applyJT() reads two input slots:  d_fout  (force at Pc_A → Beam-1)
//                                      d_foutB (force at Pc_B → Beam-2)
//
//    J shape:  (6K) × (6N1+6N2)
//    JT shape: (6N1+6N2) × (6K)
// ════════════════════════════════════════════════════════════════════════════

class BcmFixture
{
public:
    Cosserat::BeamContactMapping bcm;

    // Position Data — read by apply().
    D_In1Coord d_in1, d_in2;
    D_OutCoord d_out;    // gap mode: δ[k];  contactPoints: Pc_A = out[0][k]
    D_OutCoord d_outB;   // contactPoints mode only: Pc_B = out[1][k]

    // Velocity Data — written/read by applyJ().
    D_In1Deriv d_vel1;
    D_In2Deriv d_vel2;
    D_OutDeriv d_vout;   // gap mode or Ṗc_A
    D_OutDeriv d_voutB;  // contactPoints mode only: Ṗc_B

    // Force Data — written/read by applyJT().
    D_In1Deriv d_frc1;
    D_In2Deriv d_frc2;
    D_OutDeriv d_fout;   // gap mode or force at Pc_A
    D_OutDeriv d_foutB;  // contactPoints mode only: force at Pc_B

    int  N1{}, N2{}, K{};
    bool isContactPointsMode{false};

    /// \param mappingMode  "gap" (default) or "contactPoints"
    void setup(const In1VecCoord&        pos1,
               const In1VecCoord&        pos2,
               const std::vector<Vec2i>& sectionIds,
               const std::vector<Vec2d>& curviParams,
               bool                      algo2,
               const std::string&        mappingMode = "gap")
    {
        N1 = static_cast<int>(pos1.size());
        N2 = static_cast<int>(pos2.size());
        K  = static_cast<int>(sectionIds.size());
        isContactPointsMode = (mappingMode == "contactPoints");

        // ── Step 1: Configure BCM ────────────────────────────────────────────
        bcm.d_contactSectionIds.setValue(
            sofa::type::vector<Vec2i>(sectionIds.begin(), sectionIds.end()));
        bcm.d_curvilinearParams.setValue(
            sofa::type::vector<Vec2d>(curviParams.begin(), curviParams.end()));
        bcm.d_radius1.setValue(R1);
        bcm.d_radius2.setValue(R2);
        bcm.d_isAlgo2.setValue(algo2);
        bcm.d_mappingMode.setValue(mappingMode);

        // ── Step 2: Set input positions ──────────────────────────────────────
        d_in1.setValue(pos1);
        d_in2.setValue(pos2);

        // ── Step 3: Pre-size all Data objects ────────────────────────────────
        // apply() resizes d_out (and d_outB) itself, but the derivative Data
        // must be pre-sized because applyJ / applyJT do not resize them.
        d_out .setValue(OutVecCoord(static_cast<sofa::Size>(K)));
        d_vel1.setValue(In1VecDeriv(static_cast<sofa::Size>(N1)));
        d_vel2.setValue(In2VecDeriv(static_cast<sofa::Size>(N2)));
        d_vout.setValue(OutVecDeriv(static_cast<sofa::Size>(K)));
        d_frc1.setValue(In1VecDeriv(static_cast<sofa::Size>(N1)));
        d_frc2.setValue(In2VecDeriv(static_cast<sofa::Size>(N2)));
        d_fout.setValue(OutVecDeriv(static_cast<sofa::Size>(K)));

        if (isContactPointsMode)
        {
            d_outB .setValue(OutVecCoord(static_cast<sofa::Size>(K)));
            d_voutB.setValue(OutVecDeriv(static_cast<sofa::Size>(K)));
            d_foutB.setValue(OutVecDeriv(static_cast<sofa::Size>(K)));
        }

        // ── Step 4: Call apply() ─────────────────────────────────────────────
        // Computes output positions and populates m_jacCache.
        const sofa::core::MechanicalParams* mp =
            sofa::core::MechanicalParams::defaultInstance();

        sofa::type::vector<D_OutCoord*>       outPos = { &d_out };
        if (isContactPointsMode) outPos.push_back(&d_outB);

        sofa::type::vector<const D_In1Coord*> inPos1 = { &d_in1 };
        sofa::type::vector<const D_In1Coord*> inPos2 = { &d_in2 };

        bcm.apply(mp, outPos, inPos1, inPos2);
    }
};

// ════════════════════════════════════════════════════════════════════════════
//  Sec 5 — Jacobian probing helpers
//
//  assembleJ() and assembleJT() extract the full Jacobian matrices by
//  calling applyJ() / applyJT() with unit inputs, one DOF at a time.
//  Neither function contains any formula — they only set inputs to 1 and
//  record what the BCM outputs.
//
//  Matrix storage: row-major flat vector.  Element (r,c): mat[r*nCols + c].
//
//  gap mode:
//    J  shape: (3K)  × (6N1+6N2)
//    JT shape: (6N1+6N2) × (3K)
//
//  contactPoints mode:
//    J  shape: (6K)  × (6N1+6N2)
//              rows 0..3K-1  = Ṗc_A velocity
//              rows 3K..6K-1 = Ṗc_B velocity
//    JT shape: (6N1+6N2) × (6K)
//              cols 0..3K-1  = force at Pc_A (→ Beam-1 only)
//              cols 3K..6K-1 = force at Pc_B (→ Beam-2 only)
// ════════════════════════════════════════════════════════════════════════════

/// Zero all entries of a rigid-body VecDeriv Data (velocities or forces).
static void zeroRigidDeriv(D_In1Deriv& d, int n)
{
    d.setValue(In1VecDeriv(static_cast<sofa::Size>(n)));
}

/// Assemble J by probing applyJ() with 6N1+6N2 unit-velocity inputs.
std::vector<double> assembleJ(BcmFixture& fx)
{
    const int nRows = fx.isContactPointsMode ? 6 * fx.K : 3 * fx.K;
    const int nCols = 6 * fx.N1 + 6 * fx.N2;
    std::vector<double> J(static_cast<size_t>(nRows * nCols), 0.0);

    const sofa::core::MechanicalParams* mp =
        sofa::core::MechanicalParams::defaultInstance();

    sofa::type::vector<D_OutDeriv*>       outVels = { &fx.d_vout };
    if (fx.isContactPointsMode) outVels.push_back(&fx.d_voutB);

    sofa::type::vector<const D_In1Deriv*> in1Vels = { &fx.d_vel1 };
    sofa::type::vector<const D_In2Deriv*> in2Vels = { &fx.d_vel2 };

    // After each applyJ() call, read both output slots into column `col`.
    auto recordColumn = [&](int col)
    {
        fx.bcm.applyJ(mp, outVels, in1Vels, in2Vels);

        // Rows 0..3K-1: Pc_A velocity (or gap velocity in gap mode).
        { auto r = sofa::helper::getReadAccessor(fx.d_vout);
          for (int k = 0; k < fx.K; ++k)
              for (int d = 0; d < 3; ++d)
                  J[(3*k + d) * nCols + col] = r[k][d]; }

        // Rows 3K..6K-1: Pc_B velocity (contactPoints mode only).
        if (fx.isContactPointsMode)
        {
            auto r = sofa::helper::getReadAccessor(fx.d_voutB);
            for (int k = 0; k < fx.K; ++k)
                for (int d = 0; d < 3; ++d)
                    J[(3*fx.K + 3*k + d) * nCols + col] = r[k][d];
        }
    };

    // ── Beam-1 DOFs (cols 0 … 6N1-1) ────────────────────────────────────────
    for (int fr = 0; fr < fx.N1; ++fr) {
        for (int comp = 0; comp < 6; ++comp) {
            zeroRigidDeriv(fx.d_vel1, fx.N1);
            zeroRigidDeriv(fx.d_vel2, fx.N2);
            { auto w = sofa::helper::getWriteAccessor(fx.d_vel1);
              if (comp < 3) w[fr].getLinear()[comp]     = 1.0;
              else          w[fr].getAngular()[comp - 3] = 1.0; }
            recordColumn(6 * fr + comp);
        }
    }

    // ── Beam-2 DOFs (cols 6N1 … 6(N1+N2)-1) ─────────────────────────────────
    for (int fr = 0; fr < fx.N2; ++fr) {
        for (int comp = 0; comp < 6; ++comp) {
            zeroRigidDeriv(fx.d_vel1, fx.N1);
            zeroRigidDeriv(fx.d_vel2, fx.N2);
            { auto w = sofa::helper::getWriteAccessor(fx.d_vel2);
              if (comp < 3) w[fr].getLinear()[comp]     = 1.0;
              else          w[fr].getAngular()[comp - 3] = 1.0; }
            recordColumn(6 * fx.N1 + 6 * fr + comp);
        }
    }
    return J;
}

/// Assemble JT by probing applyJT() with unit-force inputs.
///
/// IMPORTANT: applyJT() ACCUMULATES into beam force Data via +=.
/// Beam forces MUST be zeroed before every call.
///
/// gap mode:           probes 3K columns via d_fout.
/// contactPoints mode: probes 3K columns via d_fout  (→ cols 0..3K-1)
///                     probes 3K columns via d_foutB (→ cols 3K..6K-1)
///                     The inactive slot is zeroed before each call.
std::vector<double> assembleJT(BcmFixture& fx)
{
    const int nRows = 6 * fx.N1 + 6 * fx.N2;
    const int nCols = fx.isContactPointsMode ? 6 * fx.K : 3 * fx.K;
    std::vector<double> JT(static_cast<size_t>(nRows * nCols), 0.0);

    const sofa::core::MechanicalParams* mp =
        sofa::core::MechanicalParams::defaultInstance();

    sofa::type::vector<D_In1Deriv*>       out1F = { &fx.d_frc1 };
    sofa::type::vector<D_In2Deriv*>       out2F = { &fx.d_frc2 };
    sofa::type::vector<const D_OutDeriv*> inF   = { &fx.d_fout };
    if (fx.isContactPointsMode) inF.push_back(&fx.d_foutB);

    // Harvest beam forces into column `col` of JT after applyJT().
    auto harvestColumn = [&](int col)
    {
        fx.bcm.applyJT(mp, out1F, out2F, inF);
        { auto f1 = sofa::helper::getReadAccessor(fx.d_frc1);
          for (int fr = 0; fr < fx.N1; ++fr) {
              for (int c = 0; c < 3; ++c)
                  JT[(6*fr + c)     * nCols + col] = f1[fr].getLinear()[c];
              for (int c = 0; c < 3; ++c)
                  JT[(6*fr + 3 + c) * nCols + col] = f1[fr].getAngular()[c]; }}
        { auto f2 = sofa::helper::getReadAccessor(fx.d_frc2);
          for (int fr = 0; fr < fx.N2; ++fr) {
              for (int c = 0; c < 3; ++c)
                  JT[(6*fx.N1 + 6*fr + c)     * nCols + col] = f2[fr].getLinear()[c];
              for (int c = 0; c < 3; ++c)
                  JT[(6*fx.N1 + 6*fr + 3 + c) * nCols + col] = f2[fr].getAngular()[c]; }}
    };

    // ── Columns 0 .. 3K-1: unit force in d_fout (Pc_A / gap slot) ────────────
    for (int col = 0; col < 3 * fx.K; ++col)
    {
        const int k = col / 3, d = col % 3;
        zeroRigidDeriv(fx.d_frc1, fx.N1);
        zeroRigidDeriv(fx.d_frc2, fx.N2);
        { OutVecDeriv gF(static_cast<sofa::Size>(fx.K));
          gF[k][d] = 1.0;
          fx.d_fout.setValue(gF); }
        if (fx.isContactPointsMode)
            fx.d_foutB.setValue(OutVecDeriv(static_cast<sofa::Size>(fx.K)));
        harvestColumn(col);
    }

    // ── Columns 3K .. 6K-1: unit force in d_foutB (Pc_B slot) ───────────────
    if (fx.isContactPointsMode)
    {
        for (int col = 0; col < 3 * fx.K; ++col)
        {
            const int k = col / 3, d = col % 3;
            zeroRigidDeriv(fx.d_frc1, fx.N1);
            zeroRigidDeriv(fx.d_frc2, fx.N2);
            fx.d_fout.setValue(OutVecDeriv(static_cast<sofa::Size>(fx.K)));
            { OutVecDeriv bF(static_cast<sofa::Size>(fx.K));
              bF[k][d] = 1.0;
              fx.d_foutB.setValue(bF); }
            harvestColumn(3 * fx.K + col);
        }
    }
    return JT;
}

/// Return max|J[r,c] − JT[c,r]| over all entries.
double maxSymmetryError(const std::vector<double>& J,
                         const std::vector<double>& JT,
                         int nRowsJ, int nColsJ)
{
    double maxErr = 0.0;
    for (int r = 0; r < nRowsJ; ++r)
        for (int c = 0; c < nColsJ; ++c)
            maxErr = std::max(maxErr,
                std::abs(J[r * nColsJ + c] - JT[c * nRowsJ + r]));
    return maxErr;
}

/// Return |F·(J·v) − v·(JT·F)|.
double virtualWorkError(const std::vector<double>& J,
                         const std::vector<double>& JT,
                         int nRowsJ, int nColsJ,
                         const std::vector<double>& v,
                         const std::vector<double>& F)
{
    std::vector<double> Jv(nRowsJ, 0.0);
    for (int r = 0; r < nRowsJ; ++r)
        for (int c = 0; c < nColsJ; ++c)
            Jv[r] += J[r * nColsJ + c] * v[c];
    double FJv = 0.0;
    for (int r = 0; r < nRowsJ; ++r) FJv += F[r] * Jv[r];

    std::vector<double> JTF(nColsJ, 0.0);
    for (int r = 0; r < nColsJ; ++r)
        for (int c = 0; c < nRowsJ; ++c)
            JTF[r] += JT[r * nRowsJ + c] * F[c];
    double vJTF = 0.0;
    for (int c = 0; c < nColsJ; ++c) vJTF += v[c] * JTF[c];

    return std::abs(FJv - vJTF);
}

// ════════════════════════════════════════════════════════════════════════════
//  Sec 6 — Test-case struct for the parameterised consistency tests
// ════════════════════════════════════════════════════════════════════════════

struct BcmTestCase
{
    std::string        label;
    In1VecCoord        pos1;
    In1VecCoord        pos2;
    std::vector<Vec2i> sectionIds;
    std::vector<Vec2d> curviParams;
    bool               algo2;
    std::string        mappingMode;  // "gap" or "contactPoints"
};

// ── Canonical beam geometries ─────────────────────────────────────────────────

In1VecCoord straightBeam1()
{
    return { makeRigid(0,0,0), makeRigid(1,0,0),
             makeRigid(2,0,0), makeRigid(3,0,0) };
}

In1VecCoord beam2atH()
{
    return { makeRigid(1.5,0.5,H), makeRigid(1.5,-0.5,H) };
}

In1VecCoord curvedBeam1()
{
    constexpr double R = 2.0;
    In1VecCoord v;
    for (double th : {0.0, 30.0, 60.0, 90.0}) {
        const auto [qx,qy,qz,qw] = rotZ(th);
        const double rad = th * M_PI / 180.0;
        v.push_back(makeRigid(R*std::cos(rad), R*std::sin(rad), 0.0,
                              qx, qy, qz, qw));
    }
    return v;
}

/// 16 parameterised cases: 8 geometries × 2 modes ("gap" + "contactPoints").
/// Used ONLY for consistency tests (symmetry, virtual work, sparsity,
/// no-cross-coupling).  No expected apply() output values are needed.
std::vector<BcmTestCase> makeTestCases()
{
    const auto b1s = straightBeam1();
    const auto b2s = beam2atH();
    const auto b1c = curvedBeam1();
    const auto b2c = beam2atH();

    // For each geometry, emit one gap case and one contactPoints case.
    auto both = [](const std::string&  base,
                   const In1VecCoord&  p1,
                   const In1VecCoord&  p2,
                   std::vector<Vec2i>  ids,
                   std::vector<Vec2d>  params,
                   bool                a2)
        -> std::vector<BcmTestCase>
    {
        return {
            { base + "_gap", p1, p2, ids, params, a2, "gap"           },
            { base + "_cp",  p1, p2, ids, params, a2, "contactPoints"  },
        };
    };

    std::vector<BcmTestCase> cases;
    auto add = [&](auto&& v) { for (auto& c : v) cases.push_back(c); };

    add(both("Straight_ALGO1_a050_b050",    b1s,b2s,{Vec2i{1,0}},{Vec2d{0.50,0.50}},false));
    add(both("Straight_ALGO1_a025_b075",    b1s,b2s,{Vec2i{1,0}},{Vec2d{0.25,0.75}},false));
    add(both("Straight_ALGO2_nodeC1_b075",  b1s,b2s,{Vec2i{2,0}},{Vec2d{0.00,0.75}},true ));
    add(both("Straight_ALGO2_nodeB1_b060",  b1s,b2s,{Vec2i{1,0}},{Vec2d{0.00,0.60}},true ));
    add(both("Curved_ALGO1_seg1_a050_b050", b1c,b2c,{Vec2i{1,0}},{Vec2d{0.50,0.50}},false));
    add(both("Curved_ALGO1_seg0_a040_b050", b1c,b2c,{Vec2i{0,0}},{Vec2d{0.40,0.50}},false));
    add(both("Curved_ALGO2_nodeB1_b050",    b1c,b2c,{Vec2i{1,0}},{Vec2d{0.00,0.50}},true ));
    add(both("Curved_ALGO2_nodeC1_b030",    b1c,b2c,{Vec2i{2,0}},{Vec2d{0.00,0.30}},true ));

    return cases;
}

} // anonymous namespace

// ════════════════════════════════════════════════════════════════════════════
//  Sec 7 — apply() tests — gap mode
//
//  Expected values are derived by inspection; see the file-level comment.
// ════════════════════════════════════════════════════════════════════════════

// ── ZSeparation_ALGO1 ─────────────────────────────────────────────────────────
//
// Scene
//   Beam-1: (0,0,0),(1,0,0)   Beam-2: (0.5,0.5,H),(0.5,-0.5,H)
//   α=β=0.5 → P_A=(0.5,0,0), P_B=(0.5,0,H), n̂=(0,0,1) by inspection.
//
// Expected gap = (0, 0, H−R1−R2)

TEST(Apply, ZSeparation_ALGO1)
{
    BcmFixture fx;
    fx.setup(
        { makeRigid(0,0,0), makeRigid(1,0,0) },
        { makeRigid(0.5,0.5,H), makeRigid(0.5,-0.5,H) },
        { Vec2i{0,0} }, { Vec2d{0.5,0.5} }, false);

    const auto out = sofa::helper::getReadAccessor(fx.d_out);
    ASSERT_EQ(static_cast<int>(out.size()), 1);
    EXPECT_NEAR(out[0][0], 0.0,       TOL_APPLY) << "gap.x must be 0";
    EXPECT_NEAR(out[0][1], 0.0,       TOL_APPLY) << "gap.y must be 0";
    EXPECT_NEAR(out[0][2], H-R1-R2,   TOL_APPLY) << "gap.z must be H-R1-R2";
}

// ── ZSeparation_ALGO2 ─────────────────────────────────────────────────────────
//
// Scene
//   Beam-1 node 0: (0.5,0,0)   Beam-2: (0,0,H),(1,0,H)
//   α=0 (ALGO_2), β=0.5 → P_A=(0.5,0,0), P_B=(0.5,0,H), n̂=(0,0,1).
//
// Expected gap = (0, 0, H−R1−R2)

TEST(Apply, ZSeparation_ALGO2)
{
    BcmFixture fx;
    fx.setup(
        { makeRigid(0.5,0,0), makeRigid(1.5,0,0) },
        { makeRigid(0,0,H),   makeRigid(1,0,H)   },
        { Vec2i{0,0} }, { Vec2d{0.0,0.5} }, true);

    const auto out = sofa::helper::getReadAccessor(fx.d_out);
    ASSERT_EQ(static_cast<int>(out.size()), 1);
    EXPECT_NEAR(out[0][0], 0.0,       TOL_APPLY) << "gap.x must be 0";
    EXPECT_NEAR(out[0][1], 0.0,       TOL_APPLY) << "gap.y must be 0";
    EXPECT_NEAR(out[0][2], H-R1-R2,   TOL_APPLY) << "gap.z must be H-R1-R2";
}

// ── XSeparation_ALGO1 ─────────────────────────────────────────────────────────
//
// Scene
//   Beam-1: (0,0,0),(0,0,1)   Beam-2: (D,0,0),(D,0,1)
//   α=β=0.5 → P_A=(0,0,0.5), P_B=(D,0,0.5), n̂=(1,0,0).
//
// Expected gap = (D−R1−R2, 0, 0)

TEST(Apply, XSeparation_ALGO1)
{
    BcmFixture fx;
    fx.setup(
        { makeRigid(0,0,0), makeRigid(0,0,1) },
        { makeRigid(D,0,0), makeRigid(D,0,1) },
        { Vec2i{0,0} }, { Vec2d{0.5,0.5} }, false);

    const auto out = sofa::helper::getReadAccessor(fx.d_out);
    ASSERT_EQ(static_cast<int>(out.size()), 1);
    EXPECT_NEAR(out[0][0], D-R1-R2, TOL_APPLY) << "gap.x must be D-R1-R2";
    EXPECT_NEAR(out[0][1], 0.0,     TOL_APPLY) << "gap.y must be 0";
    EXPECT_NEAR(out[0][2], 0.0,     TOL_APPLY) << "gap.z must be 0";
}

// ── ZSeparation_PenetrationIsNegative ────────────────────────────────────────
//
// Hpen=0.05 < R1+R2=0.15 → centrelines overlap → gap.z must be negative.

TEST(Apply, ZSeparation_PenetrationIsNegative)
{
    constexpr double Hpen = 0.05;
    BcmFixture fx;
    fx.setup(
        { makeRigid(0,0,0),    makeRigid(1,0,0) },
        { makeRigid(0,0,Hpen), makeRigid(1,0,Hpen) },
        { Vec2i{0,0} }, { Vec2d{0.5,0.5} }, false);

    const auto out = sofa::helper::getReadAccessor(fx.d_out);
    ASSERT_EQ(static_cast<int>(out.size()), 1);
    EXPECT_NEAR(out[0][0], 0.0, TOL_APPLY) << "gap.x must be 0";
    EXPECT_NEAR(out[0][1], 0.0, TOL_APPLY) << "gap.y must be 0";
    EXPECT_LT  (out[0][2], 0.0)            << "gap.z must be negative (penetration)";
}

// ════════════════════════════════════════════════════════════════════════════
//  Sec 8 — apply() tests — contactPoints mode
//
//  Same scenes as Sec 7.  Checks individual surface points:
//
//      Pc_A = P_A + R1·n̂   (d_out  = out[0])
//      Pc_B = P_B − R2·n̂   (d_outB = out[1])
//
//  Every test also verifies Pc_B − Pc_A == gap to confirm both modes
//  are algebraically consistent.
// ════════════════════════════════════════════════════════════════════════════

// ── ContactPoints_ZSeparation_ALGO1 ──────────────────────────────────────────
//
// Scene identical to ZSeparation_ALGO1.
// n̂=(0,0,1) by inspection.
//
// Expected
//   Pc_A = (0.5, 0, R1 )
//   Pc_B = (0.5, 0, H−R2)

TEST(Apply, ContactPoints_ZSeparation_ALGO1)
{
    BcmFixture fx;
    fx.setup(
        { makeRigid(0,0,0), makeRigid(1,0,0) },
        { makeRigid(0.5,0.5,H), makeRigid(0.5,-0.5,H) },
        { Vec2i{0,0} }, { Vec2d{0.5,0.5} }, false,
        "contactPoints");

    const auto pcA = sofa::helper::getReadAccessor(fx.d_out);   // out[0]
    const auto pcB = sofa::helper::getReadAccessor(fx.d_outB);  // out[1]

    ASSERT_EQ(static_cast<int>(pcA.size()), 1) << "K=1 Pc_A entries expected";
    ASSERT_EQ(static_cast<int>(pcB.size()), 1) << "K=1 Pc_B entries expected";

    // Pc_A = P_A + R1·n̂  where P_A=(0.5,0,0), n̂=(0,0,1).
    EXPECT_NEAR(pcA[0][0], 0.5, TOL_APPLY) << "Pc_A.x = P_A.x";
    EXPECT_NEAR(pcA[0][1], 0.0, TOL_APPLY) << "Pc_A.y = 0";
    EXPECT_NEAR(pcA[0][2], R1,  TOL_APPLY) << "Pc_A.z = R1";

    // Pc_B = P_B − R2·n̂  where P_B=(0.5,0,H), n̂=(0,0,1).
    EXPECT_NEAR(pcB[0][0], 0.5,  TOL_APPLY) << "Pc_B.x = P_B.x";
    EXPECT_NEAR(pcB[0][1], 0.0,  TOL_APPLY) << "Pc_B.y = 0";
    EXPECT_NEAR(pcB[0][2], H-R2, TOL_APPLY) << "Pc_B.z = H − R2";

    // Cross-check: Pc_B − Pc_A == gap.
    EXPECT_NEAR(pcB[0][0]-pcA[0][0], 0.0,     TOL_APPLY) << "gap.x = 0";
    EXPECT_NEAR(pcB[0][1]-pcA[0][1], 0.0,     TOL_APPLY) << "gap.y = 0";
    EXPECT_NEAR(pcB[0][2]-pcA[0][2], H-R1-R2, TOL_APPLY) << "gap.z = H−R1−R2";
}

// ── ContactPoints_ZSeparation_ALGO2 ──────────────────────────────────────────
//
// Scene identical to ZSeparation_ALGO2.
// n̂=(0,0,1).
//
// Expected
//   Pc_A = (0.5, 0, R1 )
//   Pc_B = (0.5, 0, H−R2)

TEST(Apply, ContactPoints_ZSeparation_ALGO2)
{
    BcmFixture fx;
    fx.setup(
        { makeRigid(0.5,0,0), makeRigid(1.5,0,0) },
        { makeRigid(0,0,H),   makeRigid(1,0,H)   },
        { Vec2i{0,0} }, { Vec2d{0.0,0.5} }, true,
        "contactPoints");

    const auto pcA = sofa::helper::getReadAccessor(fx.d_out);
    const auto pcB = sofa::helper::getReadAccessor(fx.d_outB);

    ASSERT_EQ(static_cast<int>(pcA.size()), 1);
    ASSERT_EQ(static_cast<int>(pcB.size()), 1);

    EXPECT_NEAR(pcA[0][0], 0.5, TOL_APPLY) << "Pc_A.x = node-0 x";
    EXPECT_NEAR(pcA[0][1], 0.0, TOL_APPLY) << "Pc_A.y = 0";
    EXPECT_NEAR(pcA[0][2], R1,  TOL_APPLY) << "Pc_A.z = R1";

    EXPECT_NEAR(pcB[0][0], 0.5,  TOL_APPLY) << "Pc_B.x";
    EXPECT_NEAR(pcB[0][1], 0.0,  TOL_APPLY) << "Pc_B.y = 0";
    EXPECT_NEAR(pcB[0][2], H-R2, TOL_APPLY) << "Pc_B.z = H − R2";

    EXPECT_NEAR(pcB[0][2]-pcA[0][2], H-R1-R2, TOL_APPLY) << "gap.z = H−R1−R2";
}

// ── ContactPoints_XSeparation_ALGO1 ──────────────────────────────────────────
//
// Scene identical to XSeparation_ALGO1.
// n̂=(1,0,0).
//
// Expected
//   Pc_A = (R1,   0, 0.5)
//   Pc_B = (D−R2, 0, 0.5)

TEST(Apply, ContactPoints_XSeparation_ALGO1)
{
    BcmFixture fx;
    fx.setup(
        { makeRigid(0,0,0), makeRigid(0,0,1) },
        { makeRigid(D,0,0), makeRigid(D,0,1) },
        { Vec2i{0,0} }, { Vec2d{0.5,0.5} }, false,
        "contactPoints");

    const auto pcA = sofa::helper::getReadAccessor(fx.d_out);
    const auto pcB = sofa::helper::getReadAccessor(fx.d_outB);

    ASSERT_EQ(static_cast<int>(pcA.size()), 1);
    ASSERT_EQ(static_cast<int>(pcB.size()), 1);

    // Pc_A = (0,0,0.5) + R1*(1,0,0)
    EXPECT_NEAR(pcA[0][0], R1,  TOL_APPLY) << "Pc_A.x = R1";
    EXPECT_NEAR(pcA[0][1], 0.0, TOL_APPLY) << "Pc_A.y = 0";
    EXPECT_NEAR(pcA[0][2], 0.5, TOL_APPLY) << "Pc_A.z = P_A.z = 0.5";

    // Pc_B = (D,0,0.5) − R2*(1,0,0)
    EXPECT_NEAR(pcB[0][0], D-R2, TOL_APPLY) << "Pc_B.x = D − R2";
    EXPECT_NEAR(pcB[0][1], 0.0,  TOL_APPLY) << "Pc_B.y = 0";
    EXPECT_NEAR(pcB[0][2], 0.5,  TOL_APPLY) << "Pc_B.z = P_B.z = 0.5";

    // Cross-check.
    EXPECT_NEAR(pcB[0][0]-pcA[0][0], D-R1-R2, TOL_APPLY) << "gap.x = D−R1−R2";
    EXPECT_NEAR(pcB[0][1]-pcA[0][1], 0.0,      TOL_APPLY) << "gap.y = 0";
    EXPECT_NEAR(pcB[0][2]-pcA[0][2], 0.0,      TOL_APPLY) << "gap.z = 0";
}

// ── ContactPoints_Radii_AreDistinct ──────────────────────────────────────────
//
// Purpose
//   Verify that R1 is applied to Pc_A and R2 to Pc_B, not the other way
//   around.  R1≠R2 means a swap shifts both points and the cross-checks
//   will catch it.
//
// Scene  (z-separation, n̂=(0,0,1))

TEST(Apply, ContactPoints_Radii_AreDistinct)
{
    BcmFixture fx;
    fx.setup(
        { makeRigid(0,0,0), makeRigid(1,0,0) },
        { makeRigid(0,0,H), makeRigid(1,0,H) },
        { Vec2i{0,0} }, { Vec2d{0.5,0.5} }, false,
        "contactPoints");

    const auto pcA = sofa::helper::getReadAccessor(fx.d_out);
    const auto pcB = sofa::helper::getReadAccessor(fx.d_outB);

    ASSERT_EQ(static_cast<int>(pcA.size()), 1);
    ASSERT_EQ(static_cast<int>(pcB.size()), 1);

    EXPECT_NEAR(pcA[0][2], R1,   TOL_APPLY) << "Pc_A.z = R1  (must NOT use R2)";
    EXPECT_NEAR(pcB[0][2], H-R2, TOL_APPLY) << "Pc_B.z = H−R2 (must NOT use R1)";
    EXPECT_NEAR(pcB[0][2]-pcA[0][2], H-R1-R2, TOL_APPLY) << "gap.z = H−R1−R2";
}

// ════════════════════════════════════════════════════════════════════════════
//  Sec 9 — applyJ translational spot-check — gap mode
//
//  Configuration: z-separation, α=β=0.5, n̂=(0,0,1).
//
//  δ̇ = Ṗc_B − Ṗc_A, rigid-body translational contribution:
//    J[gap_z, vz_beam1_frame_i  ] = −(1−α)
//    J[gap_z, vz_beam1_frame_i+1] = −α
//    J[gap_z, vz_beam2_frame_j  ] = +(1−β)
//    J[gap_z, vz_beam2_frame_j+1] = +β
//    J[gap_z, vx/vy_any_frame   ] =  0   (motion ⊥ n̂)
// ════════════════════════════════════════════════════════════════════════════

TEST(ApplyJ, TranslationalSpotCheck_ZSeparation)
{
    constexpr double alpha = 0.5, beta = 0.5;

    BcmFixture fx;
    fx.setup(
        { makeRigid(0,0,0), makeRigid(1,0,0) },
        { makeRigid(0,0,H), makeRigid(1,0,H) },
        { Vec2i{0,0} }, { Vec2d{alpha,beta} }, false);
    // N1=2, N2=2, K=1.  J shape: 3×24.

    const auto J    = assembleJ(fx);
    const int nCols = 6*fx.N1 + 6*fx.N2;  // 24
    const int B2OFF = 6*fx.N1;             // 12
    const int ROW_Z = 2;                   // gap_z row (k=0, dof z)

    // Beam-1 z-translation.
    EXPECT_NEAR(J[ROW_Z*nCols + 6*0+2], -(1.0-alpha), TOL_J)
        << "J[gap_z, vz_beam1_frame0]: weight=(1-alpha), sign=-";
    EXPECT_NEAR(J[ROW_Z*nCols + 6*1+2], -alpha, TOL_J)
        << "J[gap_z, vz_beam1_frame1]: weight=alpha, sign=-";

    // Beam-2 z-translation.
    EXPECT_NEAR(J[ROW_Z*nCols + B2OFF + 6*0+2], +(1.0-beta), TOL_J)
        << "J[gap_z, vz_beam2_frame0]: weight=(1-beta), sign=+";
    EXPECT_NEAR(J[ROW_Z*nCols + B2OFF + 6*1+2], +beta, TOL_J)
        << "J[gap_z, vz_beam2_frame1]: weight=beta, sign=+";

    // x/y motion must not affect gap_z.
    for (int fr = 0; fr < fx.N1; ++fr) {
        EXPECT_NEAR(J[ROW_Z*nCols + 6*fr+0], 0.0, TOL_J)
            << "J[gap_z, vx_beam1_" << fr << "]: ⊥ n̂, must be 0";
        EXPECT_NEAR(J[ROW_Z*nCols + 6*fr+1], 0.0, TOL_J)
            << "J[gap_z, vy_beam1_" << fr << "]: ⊥ n̂, must be 0";
    }
    for (int fr = 0; fr < fx.N2; ++fr) {
        EXPECT_NEAR(J[ROW_Z*nCols + B2OFF + 6*fr+0], 0.0, TOL_J)
            << "J[gap_z, vx_beam2_" << fr << "]: ⊥ n̂, must be 0";
        EXPECT_NEAR(J[ROW_Z*nCols + B2OFF + 6*fr+1], 0.0, TOL_J)
            << "J[gap_z, vy_beam2_" << fr << "]: ⊥ n̂, must be 0";
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  Sec 10 — applyJ translational spot-check — contactPoints mode
//
//  Configuration: z-separation, α=β=0.5, n̂=(0,0,1) frozen.
//
//  J layout: rows 0..2 = Ṗc_A,  rows 3..5 = Ṗc_B.
//
//  With n̂ frozen, Pc_A = P_A + R1·n̂ moves purely with Beam-1 (positive),
//  Pc_B = P_B − R2·n̂ moves purely with Beam-2 (positive).
//
//    J[Pc_A.z, vz_beam1_i  ] = +(1−α)   (positive — no subtraction)
//    J[Pc_A.z, vz_beam1_i+1] = +α
//    J[Pc_A.z, any_beam2   ] =  0        (no cross-coupling with frozen n̂)
//
//    J[Pc_B.z, vz_beam2_j  ] = +(1−β)
//    J[Pc_B.z, vz_beam2_j+1] = +β
//    J[Pc_B.z, any_beam1   ] =  0
// ════════════════════════════════════════════════════════════════════════════

TEST(ApplyJ, ContactPoints_TranslationalSpotCheck_ZSeparation)
{
    constexpr double alpha = 0.5, beta = 0.5;

    BcmFixture fx;
    fx.setup(
        { makeRigid(0,0,0), makeRigid(1,0,0) },
        { makeRigid(0,0,H), makeRigid(1,0,H) },
        { Vec2i{0,0} }, { Vec2d{alpha,beta} }, false,
        "contactPoints");
    // N1=2, N2=2, K=1.  J shape: 6×24.

    const auto J      = assembleJ(fx);
    const int nCols   = 6*fx.N1 + 6*fx.N2;  // 24
    const int B2OFF   = 6*fx.N1;             // 12
    const int ROW_A_Z = 2;                   // Pc_A.z row = 3*0 + 2
    const int ROW_B_Z = 3*fx.K + 2;          // Pc_B.z row = 3 + 2 = 5

    // ── Pc_A.z row ────────────────────────────────────────────────────────────
    // Positive weights; Pc_A moves with Beam-1, no sign inversion.
    EXPECT_NEAR(J[ROW_A_Z*nCols + 6*0+2], +(1.0-alpha), TOL_J)
        << "J[Pc_A.z, vz_beam1_frame0]: weight=(1-alpha), sign=+";
    EXPECT_NEAR(J[ROW_A_Z*nCols + 6*1+2], +alpha, TOL_J)
        << "J[Pc_A.z, vz_beam1_frame1]: weight=alpha, sign=+";

    // x/y Beam-1 motion ⊥ n̂.
    for (int fr = 0; fr < fx.N1; ++fr) {
        EXPECT_NEAR(J[ROW_A_Z*nCols + 6*fr+0], 0.0, TOL_J)
            << "J[Pc_A.z, vx_beam1_" << fr << "]: ⊥ n̂, must be 0";
        EXPECT_NEAR(J[ROW_A_Z*nCols + 6*fr+1], 0.0, TOL_J)
            << "J[Pc_A.z, vy_beam1_" << fr << "]: ⊥ n̂, must be 0";
    }

    // Pc_A must not depend on any Beam-2 DOF (frozen n̂).
    for (int fr = 0; fr < fx.N2; ++fr)
        for (int comp = 0; comp < 6; ++comp)
            EXPECT_NEAR(J[ROW_A_Z*nCols + B2OFF + 6*fr+comp], 0.0, TOL_J)
                << "J[Pc_A.z, beam2_fr" << fr << "_comp" << comp
                << "]: cross-coupling must be 0";

    // ── Pc_B.z row ────────────────────────────────────────────────────────────
    EXPECT_NEAR(J[ROW_B_Z*nCols + B2OFF + 6*0+2], +(1.0-beta), TOL_J)
        << "J[Pc_B.z, vz_beam2_frame0]: weight=(1-beta), sign=+";
    EXPECT_NEAR(J[ROW_B_Z*nCols + B2OFF + 6*1+2], +beta, TOL_J)
        << "J[Pc_B.z, vz_beam2_frame1]: weight=beta, sign=+";

    // x/y Beam-2 motion ⊥ n̂.
    for (int fr = 0; fr < fx.N2; ++fr) {
        EXPECT_NEAR(J[ROW_B_Z*nCols + B2OFF + 6*fr+0], 0.0, TOL_J)
            << "J[Pc_B.z, vx_beam2_" << fr << "]: ⊥ n̂, must be 0";
        EXPECT_NEAR(J[ROW_B_Z*nCols + B2OFF + 6*fr+1], 0.0, TOL_J)
            << "J[Pc_B.z, vy_beam2_" << fr << "]: ⊥ n̂, must be 0";
    }

    // Pc_B must not depend on any Beam-1 DOF.
    for (int fr = 0; fr < fx.N1; ++fr)
        for (int comp = 0; comp < 6; ++comp)
            EXPECT_NEAR(J[ROW_B_Z*nCols + 6*fr+comp], 0.0, TOL_J)
                << "J[Pc_B.z, beam1_fr" << fr << "_comp" << comp
                << "]: cross-coupling must be 0";
}

// ════════════════════════════════════════════════════════════════════════════
//  Sec 11 — Gap == Pc_B − Pc_A in Jacobian space
//
//  For identical geometry, the gap-mode Jacobian must equal the algebraic
//  difference of the contactPoints Jacobians row by row:
//
//      J_gap[row, col] == J_cp_B[3K+row, col] − J_cp_A[row, col]
//
//  This confirms both modes share the same underlying geometry and that
//  assembleJ is correct in both.
// ════════════════════════════════════════════════════════════════════════════

TEST(ApplyJ, GapEqualsContactPointsDifference)
{
    constexpr double alpha = 0.3, beta = 0.7;

    BcmFixture fxGap;
    fxGap.setup(
        { makeRigid(0,0,0), makeRigid(1,0,0) },
        { makeRigid(0,0,H), makeRigid(1,0,H) },
        { Vec2i{0,0} }, { Vec2d{alpha,beta} }, false, "gap");

    BcmFixture fxCP;
    fxCP.setup(
        { makeRigid(0,0,0), makeRigid(1,0,0) },
        { makeRigid(0,0,H), makeRigid(1,0,H) },
        { Vec2i{0,0} }, { Vec2d{alpha,beta} }, false, "contactPoints");

    const auto Jgap = assembleJ(fxGap);
    const auto Jcp  = assembleJ(fxCP);

    const int K     = fxGap.K;
    const int nCols = 6*fxGap.N1 + 6*fxGap.N2;

    for (int row = 0; row < 3*K; ++row)
        for (int col = 0; col < nCols; ++col)
        {
            const double jA   = Jcp [row         * nCols + col];
            const double jB   = Jcp [(3*K + row) * nCols + col];
            const double jgap = Jgap[row         * nCols + col];
            EXPECT_NEAR(jB - jA, jgap, TOL_J)
                << "J_cp_B[" << (3*K+row) << "," << col << "]"
                << " - J_cp_A[" << row << "," << col << "]"
                << " must equal J_gap[" << row << "," << col << "]";
        }
}

// ════════════════════════════════════════════════════════════════════════════
//  Sec 12 — Parameterised consistency tests
//
//  Run on 16 cases (8 geometries × 2 modes).  Each test checks an internal
//  algebraic identity that must hold regardless of geometry or mode.
// ════════════════════════════════════════════════════════════════════════════

class BcmConsistencyTest : public ::testing::TestWithParam<BcmTestCase>
{
protected:
    BcmFixture fx;

    void SetUp() override
    {
        const BcmTestCase& tc = GetParam();
        fx.setup(tc.pos1, tc.pos2, tc.sectionIds, tc.curviParams,
                 tc.algo2, tc.mappingMode);
    }
};

// ── Test A: Jacobian symmetry  max|J − JT.T| ≤ TOL_J ────────────────────────
//
// Purpose
//   applyJ() and applyJT() must implement adjoint (transpose) maps.
//   J[r,c] == JT[c,r] for every entry.
//
// Works in both modes because the adjoint condition is purely algebraic.
//   gap mode:           J is (3K)×(6N), JT is (6N)×(3K).
//   contactPoints mode: J is (6K)×(6N), JT is (6N)×(6K).

TEST_P(BcmConsistencyTest, JacobianSymmetry)
{
    const auto J  = assembleJ (fx);
    const auto JT = assembleJT(fx);

    const int nOut  = fx.isContactPointsMode ? 6*fx.K : 3*fx.K;
    const int nBeam = 6*fx.N1 + 6*fx.N2;

    const double err = maxSymmetryError(J, JT, nOut, nBeam);
    EXPECT_LE(err, TOL_J)
        << "max|J - JT.T| = " << err
        << "\n[case: " << GetParam().label << "]";
}

// ── Test B: Virtual-work identity  |F·(J·v) − v·(JT·F)| ≤ TOL_VW ────────────
//
// Purpose
//   Power in output space must equal power in beam space.
//   Violation means the mapping injects or extracts energy spuriously.

TEST_P(BcmConsistencyTest, VirtualWork)
{
    const auto J  = assembleJ (fx);
    const auto JT = assembleJT(fx);

    const int nOut  = fx.isContactPointsMode ? 6*fx.K : 3*fx.K;
    const int nBeam = 6*fx.N1 + 6*fx.N2;

    std::mt19937_64 rng(RNG_SEED);
    std::normal_distribution<double> dist;
    std::vector<double> v(nBeam), F(nOut);
    for (auto& x : v) x = dist(rng);
    for (auto& x : F) x = dist(rng);

    const double err = virtualWorkError(J, JT, nOut, nBeam, v, F);
    EXPECT_LE(err, TOL_VW)
        << "|F·(Jv) - v·(JT·F)| = " << err
        << "\n[case: " << GetParam().label << "]";
}

// ── Test C: ALGO_2 sparsity — exactly one Beam-1 block per contact ────────────
//
// Purpose
//   isAlgo2=true: each contact references exactly ONE Beam-1 frame (node i,
//   weight=1).  All other Beam-1 J blocks must be zero.
//
//   Checked on the Pc_A / gap rows (0..3K-1) which always carry Beam-1
//   information in both modes.

TEST_P(BcmConsistencyTest, Algo2OneBeam1Block)
{
    if (!GetParam().algo2)
        GTEST_SKIP() << "Sparsity check is for ALGO_2 only.";

    const BcmTestCase& tc = GetParam();
    const auto J     = assembleJ(fx);
    const int nCols  = 6*fx.N1 + 6*fx.N2;

    for (int k = 0; k < fx.K; ++k)
    {
        const int expected = tc.sectionIds[k][0];

        for (int fr = 0; fr < fx.N1; ++fr)
        {
            double blockMax = 0.0;
            for (int row = 3*k; row < 3*k+3; ++row)
                for (int col = 6*fr; col < 6*fr+6; ++col)
                    blockMax = std::max(blockMax, std::abs(J[row*nCols+col]));

            if (fr == expected) {
                EXPECT_GT(blockMax, TOL_J)
                    << "contact " << k << ": frame " << fr
                    << " is the contact node but its J-block is zero."
                    << "\n[case: " << tc.label << "]";
            } else {
                EXPECT_LE(blockMax, TOL_J)
                    << "contact " << k << ": frame " << fr
                    << " (expected=" << expected
                    << ") has non-zero J block (max=" << blockMax << ")."
                    << "\n[case: " << tc.label << "]";
            }
        }
    }
}

// ── Test D: contactPoints mode — no cross-beam coupling ───────────────────────
//
// Purpose
//   With n̂ frozen, Pc_A depends only on Beam-1 and Pc_B only on Beam-2.
//
//   Pc_A rows (0..3K-1)  must be zero for all Beam-2 columns.
//   Pc_B rows (3K..6K-1) must be zero for all Beam-1 columns.
//
// A non-zero entry indicates the frozen-normal linearisation was violated
// (e.g. n̂ was not actually frozen inside applyJ).

TEST_P(BcmConsistencyTest, ContactPoints_NoCrossBeamCoupling)
{
    if (!fx.isContactPointsMode)
        GTEST_SKIP() << "Cross-coupling check is for contactPoints mode only.";

    const auto J    = assembleJ(fx);
    const int nCols = 6*fx.N1 + 6*fx.N2;
    const int B2OFF = 6*fx.N1;

    // Pc_A rows must be zero in all Beam-2 columns.
    for (int k = 0; k < fx.K; ++k)
        for (int fr = 0; fr < fx.N2; ++fr) {
            double blockMax = 0.0;
            for (int row = 3*k; row < 3*k+3; ++row)
                for (int col = B2OFF + 6*fr; col < B2OFF + 6*fr+6; ++col)
                    blockMax = std::max(blockMax, std::abs(J[row*nCols+col]));
            EXPECT_LE(blockMax, TOL_J)
                << "Pc_A row, contact=" << k << ", Beam-2 frame=" << fr
                << ": cross-coupling (max=" << blockMax << ")."
                << "\n[case: " << GetParam().label << "]";
        }

    // Pc_B rows must be zero in all Beam-1 columns.
    for (int k = 0; k < fx.K; ++k)
        for (int fr = 0; fr < fx.N1; ++fr) {
            double blockMax = 0.0;
            for (int row = 3*fx.K+3*k; row < 3*fx.K+3*k+3; ++row)
                for (int col = 6*fr; col < 6*fr+6; ++col)
                    blockMax = std::max(blockMax, std::abs(J[row*nCols+col]));
            EXPECT_LE(blockMax, TOL_J)
                << "Pc_B row, contact=" << k << ", Beam-1 frame=" << fr
                << ": cross-coupling (max=" << blockMax << ")."
                << "\n[case: " << GetParam().label << "]";
        }
}

// ════════════════════════════════════════════════════════════════════════════
//  Sec 13 — Parameterised test instantiation
// ════════════════════════════════════════════════════════════════════════════

INSTANTIATE_TEST_SUITE_P(
    BcmCases, BcmConsistencyTest,
    ::testing::ValuesIn(makeTestCases()),
    [](const ::testing::TestParamInfo<BcmTestCase>& info) {
        return info.param.label;
    });

// ════════════════════════════════════════════════════════════════════════════
//  Sec 14 — main
// ════════════════════════════════════════════════════════════════════════════

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}