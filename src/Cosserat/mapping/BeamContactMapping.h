/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, development version     *
 *                (c) 2006-2024 INRIA, USTL, UJF, CNRS, MGH                   *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it    *
 * under the terms of the GNU Lesser General Public License as published by   *
 * the Free Software Foundation; either version 2.1 of the License, or (at   *
 * your option) any later version.                                             *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but        *
 * WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser    *
 * General Public License for more details.                                    *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 *******************************************************************************
 * Authors: The SOFA Team and external contributors (see Authors.txt)          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/
#pragma once
#include <Cosserat/config.h>
#include <sofa/core/Multi2Mapping.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/objectmodel/Link.h>        
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/MultiVecId.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>
#include <sofa/type/fixed_array.h>
#include <array>
#include <string>
#include <Cosserat/engine/ContactTriad.h>  

// Forward declaration – full definition included only in the .cpp.
namespace Cosserat { class SphereSweptIntersectionMethod; }

namespace sofa { namespace core { class ObjectFactory; } }

namespace Cosserat
{
    using sofa::type::Vec3d;
    using sofa::type::Vec2d;
    using sofa::type::Vec;
    using VecVec3 = sofa::type::vector<sofa::type::Vec3d>; 

    /*!
     * \class BeamContactMapping
     * @brief Multi2Mapping that maps SSIM contact-point descriptors to the Rigid3d
     *        frames of two Cosserat beams.
     *
     * Consumes:
     *   In1 – Beam-1 Rigid3d frames (N1 frames, from DiscreteCosseratMapping)
     *   In2 – Beam-2 Rigid3d frames (N2 frames, from DiscreteCosseratMapping)
     *   d_contactSectionIds  – {i,j}  from SphereSweptIntersectionMethod
     *   d_curvilinearParams  – {α,β}  from SphereSweptIntersectionMethod
     *
     * Both modes require exactly ONE connected output MechanicalObject.
     *
     * Output layout depends on d_mappingMode:
     *
     *   mode = "contactPoints":
     *     ONE output MechanicalObject of size 2K (K = number of contact pairs).
     *     Even indices  → Beam-1 surface points (Pc_A):
     *       out[0][2k]   = Pc_A[k] = P_A(α) + r₁·n̂
     *     Odd  indices  → Beam-2 surface points (Pc_B):
     *       out[0][2k+1] = Pc_B[k] = P_B(β) − r₂·n̂
     *     n̂ = (P_B − P_A) / ‖P_B − P_A‖  (fetched from SSIM, never recomputed here)
     *
     *   mode = "gap":                                                           
     *     ONE output MechanicalObject of size K.
     *     out[0][k] = Vec3(δ_n[k], δ_t1[k], δ_t2[k])  — full 3D gap in contact-local frame.
     *       δ_n  = (Pc_B − Pc_A) · n̂        (signed normal gap; < 0 = penetration)
     *       δ_t1 = (Pc_B − Pc_A) · t̂₁       (axial tangential gap)
     *       δ_t2 = (Pc_B − Pc_A) · t̂₂       (circumferential tangential gap)
     *       (to be fair, in δ_n the formula might be = (Pc_A − Pc_B) · n̂ 
     *       basically, it is the case when the two beams are nested and beam 2
     *       inner one, i.e. Pc_B is on inner tube).
     *     Contact-local frame basis (identical to SSIM d_distances convention):  
     *       n̂   – unit contact normal from SSIM (Beam-1 → Beam-2 for external;
     *              outer → inner for nested CTR)
     *       t̂₁  – normalize(τ₁ − (τ₁·n̂)·n̂), τ₁ = Beam-1 segment axial chord
     *       t̂₂  – n̂ × t̂₁  (circumferential direction)
     *
     *     applyJ  gives δ̇ = Vec3(s.Ṗ_rel·n̂, Ṗ_rel·t̂₁, Ṗ_rel·t̂₂).
     *     applyJT back-projects: d = n̂·d[0] + t̂₁·d[1] + t̂₂·d[2] before J^T.
     *
     * ALGO_1 (segment-to-segment, isAlgo2 = false):
     *   sectionIds[k] = {i,j}  →  segment i on Beam-1, segment j on Beam-2
     *   Beam-1 interpolation: frames[i]*(1−α) + frames[i+1]*α
     *   Beam-2 interpolation: frames[j]*(1−β) + frames[j+1]*β
     *
     * ALGO_2 (node-to-segment, isAlgo2 = true):
     *   SSIM guarantees that across all K contact pairs, either s1[k]≡0 for all k
     *   (Beam-1 is the node side) or s2[k]≡0 for all k (Beam-2 is the node side).
     *   BCM detects this once per apply() by scanning curvilinearParams:
     *     if ∃ k : s1[k] > 0  →  Beam-2 is node, Beam-1 is segment
     *     if ∃ k : s2[k] > 0  →  Beam-1 is node, Beam-2 is segment
     *     else (all zero)     →  degenerate; default to Beam-1 is node
     *   We scan (rather than check a single pair) because a contact landing
     *   exactly on a segment endpoint gives s=0 on the SEGMENT side too, which
     *   is not a global swap signal.
     *
     *   Node-side k:     frame[i] alone, weight = 1
     *   Segment-side k:  frames[j]*(1−γ) + frames[j+1]*γ, γ ∈ [0,1] from SSIM
     */
    class SOFA_COSSERAT_API BeamContactMapping
        : public sofa::core::Multi2Mapping<
        sofa::defaulttype::Rigid3dTypes,   ///< In1: Beam-1 frames
        sofa::defaulttype::Rigid3dTypes,   ///< In2: Beam-2 frames
        sofa::defaulttype::Vec3dTypes >     ///< Out: contact-point positions
    {
    public:
        SOFA_CLASS(BeamContactMapping,
            SOFA_TEMPLATE3(sofa::core::Multi2Mapping,
                sofa::defaulttype::Rigid3dTypes,
                sofa::defaulttype::Rigid3dTypes,
                sofa::defaulttype::Vec3dTypes));

        // ── Type aliases ─────────────────────────────────────────────────────
        using In1Types = sofa::defaulttype::Rigid3dTypes;
        using In2Types = sofa::defaulttype::Rigid3dTypes;
        using OutTypes = sofa::defaulttype::Vec3dTypes;

        using In1Coord = In1Types::Coord;
        using In2Coord = In2Types::Coord;
        using OutCoord = OutTypes::Coord;

        using In1VecCoord = In1Types::VecCoord;
        using In2VecCoord = In2Types::VecCoord;
        using OutVecCoord = OutTypes::VecCoord;

        using In1Deriv = In1Types::Deriv;
        using In2Deriv = In2Types::Deriv;
        using OutDeriv = OutTypes::Deriv;

        using In1VecDeriv = In1Types::VecDeriv;
        using In2VecDeriv = In2Types::VecDeriv;
        using OutVecDeriv = OutTypes::VecDeriv;

        using In1MatrixDeriv = In1Types::MatrixDeriv;
        using In2MatrixDeriv = In2Types::MatrixDeriv;
        using OutMatrixDeriv = OutTypes::MatrixDeriv;

        using Real = double;
        using Vec3 = sofa::type::Vec3d;
        using Vec2d = sofa::type::Vec2d;
        using Vec2i = sofa::type::Vec<2, int>;
        
        /// Output mapping mode: "contactPoints" or "gap" (default).
        ///
        /// Both modes require exactly ONE connected output MechanicalObject.
        ///
        /// "contactPoints": output size = 2K. Even indices [2k] = Pc_A (Beam-1),
        ///                  odd indices [2k+1] = Pc_B (Beam-2).
        ///                  applyJT(MatrixDeriv): passes d written by
        ///                  UnilateralLagrangianConstraint directly through J^T.
        ///
        /// "gap": output size = K.                                                
        ///   out[0][k] = Vec3(δ_n, δ_t1, δ_t2) in contact-local frame {n̂, t̂₁, t̂₂}.
        ///   δ_n  = (Pc_B − Pc_A)·n̂  (normal gap, < 0 = penetration).
        ///   δ_t1 = (Pc_B − Pc_A)·t̂₁ (axial tangential gap).
        ///   δ_t2 = (Pc_B − Pc_A)·t̂₂ (circumferential tangential gap).
        ///   (to be fair, in δ_n the formula might be = (Pc_A − Pc_B) · n̂ 
        ///   basically, it is the case when the two beams are nested and beam 2 
        ///   inner one, i.e. Pc_B is on inner tube).
        ///   Contact frame matches SSIM d_distances convention (see class doc).
        ///   applyJT(MatrixDeriv): converts Vec3(d_n, d_t1, d_t2) → d_phys = n̂·d_n + t̂₁·d_t1 + t̂₂·d_t2.
        ///   applyJT(VecDeriv):    converts Vec3(F_n, F_t1, F_t2) → F_phys = n̂·F_n + t̂₁·F_t1 + t̂₂·F_t2.
        sofa::core::objectmodel::Data<std::string>                d_mappingMode;
        
        
        /// Per-pair contact triad (n̂, t̂₁, t̂₂) written by apply().  One
        /// ContactTriad entry per contact pair k, holding the full
        /// orthonormal basis of the contact-local frame.  Downstream
        /// constraints (CPULC for unilateral + friction; any friction-aware
        /// gap-mode constraint) link here:
        ///     contactTriads = '@<bcm>.contactTriads'
        sofa::core::objectmodel::Data<VecContactTriad> d_contactTriads;
        
        ///   (Pc_B − Pc_A)·n̂ = s · δn
        /// Populated in init() from SSIM::gapSignForPublishedNormal().
        /// Downstream constraints link here:  gapSign = '@<bcm>.gapSign'
        sofa::core::objectmodel::Data<SReal>           d_gapSign;
        
        sofa::core::objectmodel::Data<VecVec3>  d_distances;       
        
        /// Mandatory link to the SphereSweptIntersectionMethod that owns contact
        /// normals. Normals are fetched via l_ssim->getContactNormal(k) and are
        /// NEVER recomputed locally in BCM.
        ///
        /// In gap mode BCM also injects the normal into applyJT (both overloads)
        /// because the new generic UnilateralLagrangianConstraint only writes
        /// Vec3(1,0,0) — it carries no geometry knowledge.
        ///
        /// In contactPoints mode the existing UnilateralLagrangianConstraint writes
        /// the physical direction; BCM passes it through J^T unchanged.
        sofa::core::objectmodel::SingleLink<
            BeamContactMapping,
            SphereSweptIntersectionMethod,
            sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK |
            sofa::core::objectmodel::BaseLink::FLAG_STOREPATH>    l_ssim;
        // ──────────────────────────────────────────────────────────────────────

        // ── Constructor / SOFA lifecycle ──────────────────────────────────────
        BeamContactMapping();
        ~BeamContactMapping() override = default;

        void init()   override;
        void reinit() override;

        // ── Multi2Mapping interface ───────────────────────────────────────────

        /// Geometric stiffness: always zero (contact normal is frozen inside apply()).
        /// Declared pure-virtual in BaseMapping; Multi2Mapping provides no default.
        void applyDJT(const sofa::core::MechanicalParams* /*mparams*/,
            sofa::core::MultiVecDerivId         /*inForce*/,
            sofa::core::ConstMultiVecDerivId    /*outForce*/) override {
        }

        /// Forward map: compute contact-point world positions and rebuild Jacobian cache.
        void apply(
            const sofa::core::MechanicalParams* mparams,
            const sofa::type::vector<sofa::core::objectmodel::Data<OutVecCoord>*>&
            dataVecOutPos,
            const sofa::type::vector<const sofa::core::objectmodel::Data<In1VecCoord>*>&
            dataVecIn1Pos,
            const sofa::type::vector<const sofa::core::objectmodel::Data<In2VecCoord>*>&
            dataVecIn2Pos) override;

        /// Tangent map (velocity propagation).
        void applyJ(
            const sofa::core::MechanicalParams* mparams,
            const sofa::type::vector<sofa::core::objectmodel::Data<OutVecDeriv>*>&
            dataVecOutVel,
            const sofa::type::vector<const sofa::core::objectmodel::Data<In1VecDeriv>*>&
            dataVecIn1Vel,
            const sofa::type::vector<const sofa::core::objectmodel::Data<In2VecDeriv>*>&
            dataVecIn2Vel) override;

        /// Transpose map – force back-propagation (VecDeriv version).
        void applyJT(
            const sofa::core::MechanicalParams* mparams,
            const sofa::type::vector<sofa::core::objectmodel::Data<In1VecDeriv>*>&
            dataVecOut1Force,
            const sofa::type::vector<sofa::core::objectmodel::Data<In2VecDeriv>*>&
            dataVecOut2Force,
            const sofa::type::vector<const sofa::core::objectmodel::Data<OutVecDeriv>*>&
            dataVecInForce) override;

        /// Transpose map – constraint Jacobian assembly (MatrixDeriv version).
        void applyJT(
            const sofa::core::ConstraintParams* cparams,
            const sofa::type::vector<sofa::core::objectmodel::Data<In1MatrixDeriv>*>&
            dataMatOut1,
            const sofa::type::vector<sofa::core::objectmodel::Data<In2MatrixDeriv>*>&
            dataMatOut2,
            const sofa::type::vector<const sofa::core::objectmodel::Data<OutMatrixDeriv>*>&
            dataMatIn) override;

    private:
        // ── Per-contact Jacobian block ────────────────────────────────────────

        /// Contribution of one parent frame to a contact-point velocity:
        ///   v_contact += weight · (v_frame + ω_frame × arm)
        struct JacBlock
        {
            int  frameIdx;  ///< global frame index in the respective input MO
            Real weight;    ///< interpolation weight: (1−α), α, (1−β), or β
            Vec3 arm;       ///< moment arm: +r₁·n̂ (Beam-1), −r₂·n̂ (Beam-2)
        };

        /// Jacobian cache for one contact pair k.
        /// ALGO_1: nBeam1Blocks = 2 (frames i and i+1).
        /// ALGO_2: nBeam1Blocks = 1 (frame i only; weight of i+1 vanishes).
        /// nBeam2Blocks = 2 for both algorithms.
        struct ContactJacEntry
        {
            sofa::type::fixed_array<JacBlock, 2> beam1Blocks;
            int nBeam1Blocks{ 0 };

            sofa::type::fixed_array<JacBlock, 2> beam2Blocks;
            int nBeam2Blocks{ 0 };

            // contact normal fetched from SSIM and frozen at apply() time.
            // Used by applyJ and applyJT to project gap velocities / forces onto n.
            // In contactPoints mode this field is unused (direction comes from ULC).
            Vec3 normal{ Vec3(0,0,1) };
            
            /// Contact-plane axial tangent t̂₁ (projected Beam-1 segment chord ⊥ n̂). 
            /// Identical to the t1_contact basis vector used by SSIM for d_distances[1].
            /// Formula: normalize(τ₁ − (τ₁·n̂)·n̂),  τ₁ = unit chord of Beam-1 segment [i, i+1].
            /// For ALGO_2: τ₁ from adjacent segment or frame local-X fallback.
            Vec3 tangent1{ Vec3(Real(1), Real(0), Real(0)) };
 
            /// Contact-plane circumferential tangent t̂₂ = n̂ × t̂₁.             
            /// Identical to the t2_contact basis vector used by SSIM for d_distances[2].
            Vec3 tangent2{ Vec3(Real(0), Real(1), Real(0)) };
            
            Real gapNormal   { Real(0) };   
            Real gapTangent1 { Real(0) };   
            Real gapTangent2 { Real(0) };  
        };

        /// Rebuilt every apply() call; consumed by applyJ / applyJT.
        sofa::type::vector<ContactJacEntry> m_jacCache;

        bool isGapMode() const { return d_mappingMode.getValue() == "gap"; }

        static constexpr Real s_eps = Real(1e-14);
    };

    void registerBeamContactMapping(sofa::core::ObjectFactory* factory);

} // namespace Cosserat