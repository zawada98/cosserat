/******************************************************************************
 * Cosserat Plugin for SOFA Framework                                         *
 *                                                                            *
 * SphereSweptIntersectionMethod.h                                            *
 *                                                                            *
 * Computes minimum distance and contact points between two Cosserat beams    *
 * modelled as sphere-swept surfaces (canal surfaces), using the centerline   *
 * Rigid3d frames produced by DiscreteCosseratMapping.                        *
 *                                                                            *
 * Theory: Lee et al., "Minimum distance between two sphere-swept surfaces",  *
 *         Computer-Aided Design, 2007.                                       *
 *                                                                            *
 * Two algorithms are provided:                                               *
 *   ALGO_1 – Segment-to-Segment (Lee-inspired)                              *
 *   ALGO_2 – Node-to-Segment    (Newton-Raphson projection, Xun-inspired)   *
 ******************************************************************************/
#pragma once

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Quat.h>
#include <sofa/type/vector.h>
#include <sofa/helper/OptionsGroup.h>

namespace Cosserat
{

    /**
     * @brief SphereSweptIntersectionMethod (SSIM)
     *
     * DataEngine component that, at every simulation step, computes the
     * minimum distance between every section pair of two Cosserat beams,
     * treating each beam as a sphere-swept (canal) surface.
     *
     * Inputs
     * ------
     *   beam1Frames   – positions of Beam 1 Rigid3d frames  (from FramesMO)
     *   beam2Frames   – positions of Beam 2 Rigid3d frames  (from FramesMO)
     *   radius1       – constant cross-section radius of Beam 1 [same unit as positions]
     *   radius2       – constant cross-section radius of Beam 2 [same unit as positions]
     *   algorithmType – "ALGO_1" (segment-to-segment) or "ALGO_2" (node-to-segment)
     *
     * Outputs (one entry per detected contact pair)
     * -------
     *   curvilinearParams   – {s1*, s2*} normalised parameter on each segment [0,1]
     *                         ALGO_2: s1 is always 0 (contact at Beam-1 node)
     *   distances           – signed gap δ = ||C1−C2|| − (r1+r2)
     *                         (negative ⇒ interpenetration)
     *   centerlinePoints1   – Pint,1  closest point on Beam 1 centreline
     *   centerlinePoints2   – Pint,2  closest point on Beam 2 centreline
     *   surfacePoints1      – Psurf,1 as Rigid3d:
     *                           position    = Pint,1 + r1 * n̂
     *                           orientation = contact frame (SLERP of bracketing Beam-1 frames)
     *   surfacePoints2      – Psurf,2 as Rigid3d:
     *                           position    = Pint,2 - r2 * n̂
     *                           orientation = contact frame (SLERP of bracketing Beam-2 frames)
     *   contactSectionIds   – {i, j} section/node index pair per contact
     *                         ALGO_1: i = Beam-1 segment index, j = Beam-2 segment index
     *                         ALGO_2: i = Beam-1 node  index,  j = Beam-2 segment index
     *
     * Contact-frame convention
     * ------------------------
     *   Z = n̂  (contact normal, pointing from Pint2 → Pint1)
     *   X = beam tangent projected onto the contact plane
     *   Y = Z × X  (bi-normal, right-hand rule)
     *
     * Usage in a SOFA scene (Python)
     * --------------------------------
     *   node.addObject('SphereSweptIntersectionMethod',
     *                  name          = 'ssim',
     *                  beam1Frames   = beam1.FramesMO.getLinkPath(),
     *                  beam2Frames   = beam2.FramesMO.getLinkPath(),
     *                  radius1       = 0.15,
     *                  radius2       = 0.15,
     *                  algorithmType = 'ALGO_1')
     */
    namespace
    {
        using sofa::Data;
    }

    class SphereSweptIntersectionMethod : public sofa::core::DataEngine
    {
    public:
        SOFA_CLASS(SphereSweptIntersectionMethod, sofa::core::DataEngine);

        // ── Type aliases ─────────────────────────────────────────────────────────
        using Rigid3dTypes = sofa::defaulttype::Rigid3dTypes;
        using RigidCoord = Rigid3dTypes::Coord;       ///< position + quaternion
        using VecRigidCoord = Rigid3dTypes::VecCoord;
        using Real = sofa::type::Vec3d::value_type;
        using Vec3 = sofa::type::Vec3d;
        using VecVec3 = sofa::type::vector<Vec3>;
        using VecReal = sofa::type::vector<Real>;
        using Vec2i = sofa::type::Vec<2, int>;
        using VecVec2i = sofa::type::vector<Vec2i>;
        using Vec6 = sofa::type::Vec<6, Real>;
        using VecVec6 = sofa::type::vector<Vec6>;

        // ── Inputs ───────────────────────────────────────────────────────────────
        /// Rigid3d frames of Beam 1 (output of DiscreteCosseratMapping)
        Data<VecRigidCoord> d_beam1Frames;
        /// Rigid3d frames of Beam 2 (output of DiscreteCosseratMapping)
        Data<VecRigidCoord> d_beam2Frames;
        /// Cross-section radius of Beam 1 (metres / same unit as scene)
        Data<Real>          d_radius1;
        /// Cross-section radius of Beam 2
        Data<Real>          d_radius2;
        /// "ALGO_1" = segment-to-segment (Lee),  "ALGO_2" = node-to-segment (Xun)
        Data<sofa::helper::OptionsGroup> d_algorithmType;
        /// Maximum number of Newton–Raphson iterations (ALGO_2 only)
        Data<int>           d_maxNRIterations;
        /// Newton–Raphson convergence tolerance
        Data<Real>          d_nrTolerance;

        // ── Outputs ──────────────────────────────────────────────────────────────
        /// Normalised curvilinear parameters {s1*, s2*} for each contact pair.
        /// ALGO_2: s1 = 0 always (contact located exactly at Beam-1 node i).
        Data<sofa::type::vector<sofa::type::Vec2d>> d_curvilinearParams;

        /// Gap vector {δn, δt1, δt2} per contact pair.               
        /// δn = ||Pint1−Pint2|| − (r1+r2)  (negative ⇒ interpenetration).
        /// Tangential components δt1, δt2 are reserved and set to zero.
        Data<VecVec3>       d_distances;

        /// Closest point on Beam 1 centreline per contact pair
        Data<VecVec3>       d_centerlinePoints1;
        /// Closest point on Beam 2 centreline per contact pair
        Data<VecVec3>       d_centerlinePoints2;

        //   position    = Psurf = Pint ± r * n̂
        //   orientation = contact frame built from SLERP of the two bracketing beam
        //                 frames at parameter s (see contact-frame convention above).
        /// Surface contact point on Beam 1 as Vec3d
        Data<VecVec3> d_surfacePoints1;
        /// Surface contact point on Beam 2 as Vec3d
        Data<VecVec3> d_surfacePoints2;

        /// Indices {i, j} of the beam sections (ALGO_1) or node/section (ALGO_2)
        /// forming each contact pair
        Data<VecVec2i>      d_contactSectionIds;

        // ── Constructor / SOFA lifecycle ─────────────────────────────────────────
        SphereSweptIntersectionMethod();
        ~SphereSweptIntersectionMethod() override = default;

        void init()     override;
        void reinit()   override;
        void doUpdate() override;

    protected:
        // ── Internal geometry helpers ─────────────────────────────────────────────

        /**
         * @brief ALGO_1 – Segment-to-Segment minimum distance.
         *
         * Implements the Ericson algorithm (Real-Time Collision Detection, 2005)
         * to find the pair (s1*, s2*) in [0,1]² minimising ||L1(s1)−L2(s2)||.
         *
         * @param[in]  p0,p1   endpoints of segment 1
         * @param[in]  q0,q1   endpoints of segment 2
         * @param[out] s1      parameter on segment 1
         * @param[out] s2      parameter on segment 2
         * @param[out] cp1     closest point on segment 1
         * @param[out] cp2     closest point on segment 2
         * @return true on success
         */
        static bool segmentToSegment(const Vec3& p0, const Vec3& p1,
            const Vec3& q0, const Vec3& q1,
            Real& s1, Real& s2,
            Vec3& cp1, Vec3& cp2);

        /**
         * @brief ALGO_2 – Node-to-Segment minimum distance via Newton–Raphson.
         *
         * For node P on Beam 1, find s2* on the segment Q(s)=q0+s*(q1−q0) such
         * that d/ds || P − Q(s) ||² = 0.  For a linear segment this reduces to a
         * single-step projection; NR is provided for extensibility to higher-order
         * centrelines.
         *
         * @param[in]  node    point on Beam 1 (frame position)
         * @param[in]  q0,q1   endpoints of the candidate segment on Beam 2
         * @param[out] s2      minimising parameter ∈ [0,1]
         * @param[out] cp      closest point on segment
         * @param[in]  maxIter maximum NR iterations
         * @param[in]  tol     convergence tolerance
         */
        static bool nodeToSegmentNR(const Vec3& node,
            const Vec3& q0, const Vec3& q1,
            Real& s2, Vec3& cp,
            int maxIter = 20, Real tol = 1e-12);

        /**
         * @brief Map centreline contact points to Rigid3d surface contact frames.
         *         *
         * Position:
         *   psurf1 = pint1 + r1 * n̂,   psurf2 = pint2 - r2 * n̂
         *   where n̂ = (pint2 − pint1) / ||pint2 − pint1||
         *
         * Orientation (contact frame, same convention for both surface points):
         *   Z = n̂  (contact normal)
         *   X = SLERP-interpolated beam tangent projected onto the contact plane
         *   Y = Z × X  (right-hand bi-normal)
         *   Quaternion computed via SLERP of the two bracketing beam frames at s.
         *
         * Fallback when ||pint1−pint2|| < ε:
         *   n̂ = frameA local Y-axis (coincident centrelines).
         *
         * @param pint1,pint2   centreline contact points
         * @param s1,s2         local parameters in [0,1] within each segment
         * @param r1,r2         beam radii
         * @param psurf1,psurf2 output surface contact 
         */
        static void mapToSurface(const Vec3& pint1,
            const Vec3& pint2,
            const RigidCoord& frameA,
            Real r1, Real r2,
            Vec3& psurf1,
            Vec3& psurf2);

        // ── Broad-phase candidate selection ──────────────────────────────────────
        /**
         * @brief For section i on Beam 1 return the indices of candidate sections
         *        on Beam 2 (currently returns all segments – no false negatives).
         *        Replace with a BVH / spatial hash for large N.
         */
        static sofa::type::vector<int>
            candidateSections(int i,
                const VecRigidCoord& frames1,
                const VecRigidCoord& frames2);

    private:
        static constexpr Real s_eps = Real(1e-14); ///< numerical zero
    };

} // namespace Cosserat