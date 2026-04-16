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
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/MultiVecId.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>
#include <sofa/type/fixed_array.h>
#include <array>
#include <string>

namespace sofa { namespace core { class ObjectFactory; } }

namespace Cosserat
{
    namespace
    {
        using sofa::type::Vec3d;
        using sofa::type::Vec2d;
        using sofa::type::Vec;
    } // anonymous namespace

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
     *     n̂ = (P_B − P_A) / ‖P_B − P_A‖
     *
     *   mode = "gap":
     *     ONE output MechanicalObject of size K.
     *     out[0][k] = δ[k] = Pc_B[k] − Pc_A[k] = (d − r₁ − r₂)·n̂
     *     δ > 0 : separation;  δ < 0 : penetration.
     *
     * ALGO_1 (segment-to-segment, isAlgo2 = false):
     *   sectionIds[k] = {i,j}  →  segment i on Beam-1, segment j on Beam-2
     *   Beam-1 interpolation: frames[i]*(1−α) + frames[i+1]*α
     *   Beam-2 interpolation: frames[j]*(1−β) + frames[j+1]*β
     *
     * ALGO_2 (node-to-segment, isAlgo2 = true):
     *   sectionIds[k] = {i,j}  →  node i on Beam-1, segment j on Beam-2
     *   Beam-1: frame[i] alone, weight = 1  (α = 0 always)
     *   Beam-2: frames[j]*(1−β) + frames[j+1]*β
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

        // ── Data fields (SSIM inputs) ─────────────────────────────────────────
        /// Per-contact {i,j} index pair: Beam-1 segment (ALGO_1) or node (ALGO_2) index i,
        /// Beam-2 segment index j.  Produced by SphereSweptIntersectionMethod.
        sofa::core::objectmodel::Data<sofa::type::vector<Vec2i>>  d_contactSectionIds;

        /// Per-contact {α,β} curvilinear parameters on [0,1].
        /// ALGO_2: α = 0 always (contact at Beam-1 node i).
        sofa::core::objectmodel::Data<sofa::type::vector<Vec2d>>  d_curvilinearParams;

        /// Cross-section radius of Beam-1 (same length unit as frame positions).
        sofa::core::objectmodel::Data<Real>                       d_radius1;

        /// Cross-section radius of Beam-2.
        sofa::core::objectmodel::Data<Real>                       d_radius2;

        /// True when SSIM runs ALGO_2 (node-to-segment); false for ALGO_1 (segment-to-segment).
        sofa::core::objectmodel::Data<bool>                       d_isAlgo2;

        /// Output mapping mode: "contactPoints" or "gap" (default).
        ///
        /// Both modes require exactly ONE connected output MechanicalObject.
        ///
        /// "contactPoints": output size = 2K. Even indices [2k] = Pc_A (Beam-1),
        ///                  odd indices [2k+1] = Pc_B (Beam-2).
        /// "gap"          : output size = K.  Index [k] = Pc_B[k] - Pc_A[k].
        sofa::core::objectmodel::Data<std::string>                d_mappingMode;

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
        };

        /// Rebuilt every apply() call; consumed by applyJ / applyJT.
        sofa::type::vector<ContactJacEntry> m_jacCache;

        bool isGapMode() const { return d_mappingMode.getValue() == "gap"; }

        static constexpr Real s_eps = Real(1e-14);
    };

    void registerBeamContactMapping(sofa::core::ObjectFactory* factory);

} // namespace Cosserat