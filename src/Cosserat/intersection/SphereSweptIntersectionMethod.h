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
 *         Computer-Aided Design, 2007.                                       *                                                                       *
 ******************************************************************************/
#pragma once

#include <Cosserat/config.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Quat.h>
#include <sofa/type/vector.h>
#include <sofa/helper/OptionsGroup.h>
#include <map>
#include <utility>

namespace Cosserat
{
	
	/**
	* @brief SphereSweptIntersectionMethod (SSIM)
	*
	* BaseObject service used by BeamContactMapping to compute the minimum centreline
	* distance and contact geometry for every section pair of two Cosserat beams.
	* Consumers explicitly request ContactEvaluation snapshots through evaluateContacts();
	* SSIM does not publish per-contact SOFA Data.
	* ─── Data inputs ────────────────────────────────────────────────────────
	*   beam1Frames, beam2Frames                     Rigid3d frames (FramesMO output
	*                                                of DiscreteCosseratMapping)
	*   beam1Velocities, beam2Velocities             FramesMO velocities. MANDATORY
	*                                                (see validateParameters()).
	*                                                Link via:
	*                                                  beam1Velocities =
	*                                                    beam1_MO.getLinkPath()
	*                                                    + '.velocity'
	*   radius1,  radius2                            Outer cross-section radii (>0)
	*   innerRadius1, innerRadius2                   Bore radii. 0 = solid beam.
	*   contactConfiguration                         "external" | "nested"
	*   defaultNormal                                Last-resort normal for
	*                                                coincident + parallel beams.
	*                                                MANDATORY in nested mode.
	*   broadPhaseMarginFactor                       Bounding-sphere inflation
	*                                                (default 1.5).
	*   cachedNormalMaxAxialProjection               Reject cached normal if
	*                                                |n·tangent| exceeds this
	*                                                (default sin 10°).
	*   cachedNormalMaxFrameRotation                 Reject cached normal if the
	*                                                slerp'd contact-point frame
	*                                                rotated more than this
	*                                                (default π/18 = 10°).
	*
	* evaluateContacts() returns distances, centreline/surface points,
	* section ids, contact triads, and curvilinear parameters for all contact pairs.
	* ─── Gap equations ──────────────────────────────────────────────────────
	*   External (both beams solid, side-by-side):
	*       δn = dist − (r1 + r2)
	*   Nested (one beam hollow, one solid, coaxial):
	*       δn = ri_outer − r_outer_inner − dist
	*   In both, positive δn means clearance, negative means penetration.
	*   validateParameters() verifies ri_outer > r_outer_inner at load time.
	*
	* ─── Contact normal convention ──────────────────────────────────────────
	*   computeContactNormal() internally produces n̂_raw = (Pb − Pa) / ‖·‖,
	*   i.e. Beam-1 → Beam-2. evaluateContacts() then publishes n̂_out:
	*       external:             n̂_out = n̂_raw
	*       nested, Beam1 outer:  n̂_out = n̂_raw           (outer → inner)
	*       nested, Beam1 inner:  n̂_out = −n̂_raw          (outer → inner)
	*   So the published normal always points outer → inner in nested mode,
	*   independent of scene-level beam labelling.
	*
	* ─── Gap sign (relation between δn and n̂_out) ──────────────────────────
	*   (Pc_B − Pc_A) · n̂_out = gapSignForPublishedNormal() · δn
	*       external:             +1
	*       nested, Beam1 outer:  −1
	*       nested, Beam1 inner:  +1
	*   Callers that need forces/velocities consistent with δn must multiply
	*   (vcB − vcA) · n̂_out by this sign (BeamContactMapping does this in
	*   applyJ / applyJT for the normal component only — tangentials carry no
	*   sign convention).
	*
	* ─── Tangential slip (δt1, δt2) ─────────────────────────────────────────
	*   t̂1 = normalize(τ1 − (τ1·n̂_out)·n̂_out),  τ1 = Beam-1 segment chord
	*   t̂2 = n̂_out × t̂1
	*   δt1 = (vc_Pb − vc_Pa) · t̂1 · dt
	*   δt2 = (vc_Pb − vc_Pa) · t̂2 · dt
	*   where vc_Pa / vc_Pb include both linear and angular contributions of the
	*   bounding frames (v_centre + ω × arm, blended by the curvilinear
	*   parameter).
	*
	* ─── Zero-normal fallback (fired inside computeContactNormal) ───────────
	*   When ‖Pb − Pa‖ < ε the normal cannot be derived from the centrelines.
	*   The fallback chain is, in order:
	*     1. Reuse m_lastValidNormal[(i,j)] IF it still passes both:
	*        (a) perpendicularity — |n·seg1tan|,|n·seg2tan| < axial tolerance
	*        (b) smoothness       — slerp'd frame rotation since cache time
	*                               < rotation tolerance
	*     2. normalize(seg1tan × seg2tan)  (fails if tangents are parallel)
	*     3. defaultNormal rotated by the current slerp'd Beam-1 frame
	*        (requires defaultNormal ≠ 0 — enforced in nested mode)
	*   Nested CTR beams are coaxial by construction, so the fallback is the
	*   hot path there; defaultNormal is mandatory.
	*
	*
	* ─── Python scene usage ─────────────────────────────────────────────────
	*   node.addObject('SphereSweptIntersectionMethod',
	*       name                   = 'ssim',
	*       beam1Frames            = beam1.FramesMO.getLinkPath(),
	*       beam2Frames            = beam2.FramesMO.getLinkPath(),
	*       beam1Velocities        = beam1.FramesMO.getLinkPath() + '.velocity',
	*       beam2Velocities        = beam2.FramesMO.getLinkPath() + '.velocity',
	*       radius1                = 0.15,
	*       radius2                = 0.10,
	*       innerRadius1           = 0.12,     # hollow outer tube
	*       innerRadius2           = 0.0,      # solid inner tube
	*       contactConfiguration   = 'nested',
	*       defaultNormal          = '0 1 0',  # REQUIRED for nested
    *       broadPhaseMarginFactor = 1.5)
    */	
		
    class SOFA_COSSERAT_API SphereSweptIntersectionMethod : public sofa::core::objectmodel::BaseObject
    {
    public:
        SOFA_CLASS(SphereSweptIntersectionMethod, sofa::core::objectmodel::BaseObject);

        // ── Type aliases ─────────────────────────────────────────────────────────
        using Rigid3dTypes  = sofa::defaulttype::Rigid3dTypes;
        using RigidCoord    = Rigid3dTypes::Coord;
        using RigidDeriv    = Rigid3dTypes::Deriv;               
        using VecRigidDeriv = Rigid3dTypes::VecDeriv;   
        using VecRigidCoord = Rigid3dTypes::VecCoord;
        using Real          = sofa::type::Vec3d::value_type;
        using Vec3          = sofa::type::Vec3d;
        using VecVec3       = sofa::type::vector<Vec3>;
        using VecReal       = sofa::type::vector<Real>;
        using Vec2i         = sofa::type::Vec<2, int>;
        using VecVec2i      = sofa::type::vector<Vec2i>;
		using Vec2d         = sofa::type::Vec2d;
        using VecVec2d      = sofa::type::vector<Vec2d>;

        struct ContactEvaluation
        {
            VecVec3  distances;
            VecVec3  centerlinePoints1;
            VecVec3  centerlinePoints2;
            VecVec3  surfacePoints1;
            VecVec3  surfacePoints2;
            VecVec2i contactSectionIds;
            VecVec3  contactNormals;
            VecVec3  contactTangents1;
            VecVec3  contactTangents2;
            VecVec2d curvilinearParams;
        };

        // ── Inputs ───────────────────────────────────────────────────────────────

        /// Rigid3d frames of Beam 1 (output of DiscreteCosseratMapping / FramesMO)
        sofa::Data<VecRigidCoord> d_beam1Frames;
        /// Rigid3d frames of Beam 2
        sofa::Data<VecRigidCoord> d_beam2Frames;
        
        /// Rigid3d frame velocities of Beam 1 (FramesMO.velocity).
        /// Size must match beam1Frames; init() marks the component Invalid otherwise.
        /// Link in Python: beam1Velocities = beam1_MO.getLinkPath() + '.velocity'
        sofa::Data<VecRigidDeriv> d_beam1Velocities;                   
        /// Rigid3d frame velocities of Beam 2. Same as above.
        sofa::Data<VecRigidDeriv> d_beam2Velocities;  
    	
        /// Outer cross-section radius of Beam 1
        sofa::Data<Real>          d_radius1;
        /// Outer cross-section radius of Beam 2
        sofa::Data<Real>          d_radius2;

        /// Inner radius of Beam 1. 0 = solid beam (external contact).
        /// >0 = hollow tube; Beam 1 is the CTR outer tube if radius1 > radius2.
        sofa::Data<Real>          d_innerRadius1;
        /// Inner radius of Beam 2. 0 = solid beam (external contact).
        /// >0 = hollow tube; Beam 2 is the CTR outer tube if radius2 > radius1.
        sofa::Data<Real>          d_innerRadius2;

        /// Contact geometry configuration.
        ///   "external" – beams are always side-by-side.
        ///                Gap = dist − (r1 + r2).
        ///   "nested"   – beams are always coaxial (CTR). The beam with the
        ///                larger outer radius is the outer tube.
        ///                Gap = ri_outer − r_outer_inner − dist.
        ///                A one-time init check verifies ri_outer > r_outer_inner.
        sofa::Data<sofa::helper::OptionsGroup> d_contactConfiguration;

        /// Multiplier applied to each beam's contact-relevant radius to build the
        /// bounding sphere used in broad-phase culling.
        ///   External: R_i = halfLength_i + factor × r_outer_i
        ///   Nested (outer tube): R_i = halfLength_i + factor × ri_outer_i  (bore)
        ///   Nested (inner tube): R_i = halfLength_i + factor × r_outer_i
        /// Default 1.5 catches all segments within 1.5 tube-radii of touching.
        sofa::Data<Real>          d_broadPhaseMarginFactor;


        /// Default contact normal used when centrelines are coincident AND no previous
        /// valid normal exists AND tangent cross-product is zero (parallel beams).
        /// Must be set by the user. No hardcoded fallback is applied;
        /// if this vector is zero, an error is emitted and the pair is skipped.
        /// REQUIRED for nested CTR scenes (centrelines nearly coincident by design).
        sofa::Data<Vec3> d_defaultNormal; 
        
        /// Tolerance on |n_cached · t̂| for the cached-normal perpendicularity test.
        /// The tangents used are the raw SEGMENT CHORD tangents (normalise(P_{i+1}-P_i)),
        /// not the contact-plane tangents t̂₁/t̂₂ — the latter are constructed
        /// orthogonal to n̂ and would make the test tautological.
        /// A cached normal is rejected if |n·seg1Tangent| or |n·seg2Tangent|
        /// exceeds this value. Default 0.17 ≈ sin(10°).
        sofa::Data<Real> d_cachedNormalMaxAxialProjection;       
        
        // ── Constructor / SOFA lifecycle ─────────────────────────────────────────
        SphereSweptIntersectionMethod();
        ~SphereSweptIntersectionMethod() override = default;

        void init()     override;
        void reinit()   override;

        // ── Public API for external C++ components (e.g. BeamContactMapping) ─────
        ContactEvaluation evaluateContacts(const sofa::Data<VecRigidCoord>& frames1Data,
                                           const sofa::Data<VecRigidCoord>& frames2Data,
                                           const sofa::Data<VecRigidDeriv>& vels1Data,
                                           const sofa::Data<VecRigidDeriv>& vels2Data);
        int getEvaluationParametersCounter() const;
        /**
         * @brief Computes and returns the contact normal n̂ for a given pair of
         *        centreline contact points, updating m_lastValidNormal[(i,j)] as a
         *        side effect.
         *
         * This is the canonical normal-computation entry point. It is called
         * internally by evaluateContacts() and can also be called from BeamContactMapping
         * (via a pointer to this SSIM object) if the mapping needs to recompute
         * the normal outside of the last published evaluation.
         *
         * Returns the RAW Beam-1 → Beam-2 direction (before any nested sign flip).
         * The sign flip to enforce outer→inner convention is applied by evaluateContacts()
         * on the output, not inside this function.
         *
         * @param pint1       Closest point on Beam 1 centreline
         * @param pint2       Closest point on Beam 2 centreline
         * @param i           Beam 1 section index (for cache key)
         * @param j           Beam 2 section index (for cache key)
         * @param seg1Tangent Unit tangent of Beam-1 segment i  (= normalise(P_{i+1}−P_i))
         * @param seg2Tangent Unit tangent of Beam-2 segment j  (= normalise(Q_{j+1}−Q_j))
         * @param frameA      Beam-1 frame at index i   (for SLERP fallback)
         * @param frameB      Beam-1 frame at index i+1 (for SLERP fallback)
         * @param s1          Interpolation parameter in [0,1] for SLERP fallback
         * @return            Unit contact normal n̂ (raw Beam-1 → Beam-2 direction)
         */
        Vec3 computeContactNormal(const Vec3&       pint1,
                                  const Vec3&       pint2,
                                  int               i,
                                  int               j,
                                  const Vec3&       seg1Tangent,
                                  const Vec3&       seg2Tangent,
                                  const RigidCoord& frameA,
                                  const RigidCoord& frameB,
                                  Real              s1);

        /**
         * @brief Returns the sign factor that converts SSIM's published
         *        signed-clearance gap (d_distances[k][0] = bestGap, positive = clear)
         *        into the convention  δ_n = (Pc_B − Pc_A) · nHat_out.
         *
         * Needed by downstream mappings (e.g. BeamContactMapping in gap mode)
         * that pair the published gap with the published normal in applyJ /
         * applyJT. The identity (Pc_B − Pc_A)·nHat_out = signFactor · bestGap
         * holds in all configurations:
         *
         *   External:            +1
         *   Nested, Beam1=outer: −1   (nHat_out = Beam-1→Beam-2, geometry flips)
         *   Nested, Beam1=inner: +1
         *
         * This method exists solely to keep BCM decoupled from SSIM's internal
         * nesting classification: BCM multiplies d_distances[k][0] by this sign
         * and gets a scalar consistent with its own applyJ output.
         *
         * @return +1 or −1.
         */
        Real gapSignForPublishedNormal() const; 

    protected:
        // ── Internal geometry helpers ─────────────────────────────────────────────
    	
    	/**                                                                          
		 * @brief Nested-mode segment-pair contact metric.                           
		 *                                                                           
		 * For a candidate T2 segment [q0,q1] against an outer-loop T1 segment       
		 * [p0,p1], finds the parameter range on T1 that T2 actually shadows         
		 * axially, and returns the point of MAXIMUM radial separation over that     
		 * sub-range. Replaces standard min-centreline-distance selection for        
		 * nested geometry, where worst penetration corresponds to LARGEST radial    
		 * offset (gap = ri_outer − r_outer_inner − dist).                           
		 *                                                                           
		 * Theory:                                                                   
		 *   At each s_1 ∈ [s_start, s_end] on T1, define the corresponding T2      
		 *   point Q(s_1) by enforcing (Q − P(s_1)) ⊥ d1. Then                       
		 *     s_2(s_1) = (s_1·|d1|² − (Q0−P0)·d1) / (d2·d1)                         
		 *   Q − P stays in the plane perpendicular to d1 and is linear in s_1,      
		 *   so r²(s_1) is convex quadratic → max is at an endpoint of the overlap.  
		 *                                                                           
		 * Axial-overlap filter: T2 segments whose s_1 projection lies entirely      
		 * outside [0,1] are physically in free space relative to this T1 segment    
		 * and are rejected (return false). Without this filter, the nested gap      
		 * formula produces spurious large penetrations from arc tips that hang      
		 * outside the bore axially.                                                 
		 *                                                                           
		 * @param p0,p1     Endpoints of outer-loop (T1) segment.                    
		 * @param q0,q1     Endpoints of candidate (T2) segment.                     
		 * @param s1_out    [out] Parameter on T1 at max-radial location ∈ [0,1].   
		 * @param s2_out    [out] Parameter on T2 at corresponding point ∈ [0,1].   
		 * @param cp1,cp2   [out] Contact points on T1 and T2 (P, Q).                
		 * @param radial    [out] ‖cp2 − cp1‖ at the max-radial location.            
		 * @return          true if axial overlap is non-empty; false → skip pair.   
		 */                                                                          
		static bool axialOverlapMaxRadial(                                           
		    const Vec3& p0, const Vec3& p1,                                          
		    const Vec3& q0, const Vec3& q1,                                          
		    Real& s1_out, Real& s2_out,                                              
		    Vec3& cp1, Vec3& cp2,                                                    
		    Real& radial);  

        /**         *
         * @param p0,p1   endpoints of segment 1
         * @param q0,q1   endpoints of segment 2
         * @param s1      [out] parameter on segment 1
         * @param s2      [out] parameter on segment 2
         * @param cp1     [out] closest point on segment 1
         * @param cp2     [out] closest point on segment 2
         * @return true always (degenerate cases handled internally)
         */
        static bool segmentToSegment(const Vec3& p0, const Vec3& p1,
                                     const Vec3& q0, const Vec3& q1,
                                     Real& s1, Real& s2,
                                     Vec3& cp1, Vec3& cp2);
    	
    	/**
		*
		* For a LINEAR segment Q(s) = q0 + s*(q1−q0), the minimiser of
		* ||P − Q(s)||² is given in closed form by:
		*   s* = clamp( (P−q0)·(q1−q0) / ||q1−q0||² , 0, 1 )
		*
		* @param node   point on Beam 1
		* @param q0,q1  endpoints of the candidate Beam-2 segment
		* @param s2     [out] minimising parameter ∈ [0,1]
		* @param cp     [out] closest point on the segment
		* @return true always
		*/
    	static bool nodeToSegment(const Vec3& node,
								  const Vec3& q0, const Vec3& q1,
								  Real& s2, Vec3& cp);


        /**
         * @brief Maps centreline contact points to physical surface contact points.
         *
         * The callers are responsible for passing radii appropriate to the contact
         * configuration:
         *
         *   External contact:
         *     r_surf1 = r_outer_beam1, r_surf2 = r_outer_beam2
         *   Nested contact, Beam-1 = outer tube:
         *     r_surf1 = ri_beam1 (bore),  r_surf2 = r_outer_beam2
         *   Nested contact, Beam-1 = inner tube:
         *     r_surf1 = r_outer_beam1,    r_surf2 = ri_beam2 (bore)
         *
         * In all cases the formula uses the raw Beam-1→Beam-2 normal nHat:
         *   psurf1 = pint1 +/- r_surf1 * nHat
         *   psurf2 = pint2 +/− r_surf2 * nHat
         *
         * @param pint1,pint2    centreline contact points
         * @param nHat           pre-computed unit contact normal (raw Beam-1 → Beam-2)
         * @param r_surf1        contact-relevant radius for Beam 1
         * @param r_surf2        contact-relevant radius for Beam 2
         * @param psurf1,psurf2  [out] physical surface contact positions
         */
        static void mapToSurface(const Vec3& pint1,
                                 const Vec3& pint2,
                                 const Vec3& nHat,
                                 Real  r_surf1, Real r_surf2,
                                 Vec3& psurf1, Vec3& psurf2);

        /**
         * @brief Broad-phase candidate selection 
         *
         * Tests whether the bounding sphere of the query segment 
         * overlaps the bounding sphere of each candidate segment j.
         *
         *   R_i = queryHalfLength + factor × r_bp_query
         *   R_j = halfLen_j       + factor × r_bp_candidate
         *   Segment j is a candidate iff  ||queryMid − mid_j|| ≤ R_i + R_j
         *
         * The callers supply contact-relevant broad-phase radii:
         *   External:  r_bp_query = r_outer_finer,  r_bp_candidate = r_outer_coarser
         *   Nested (finer=outer): r_bp_query = ri_finer (bore), r_bp_candidate = r_outer_coarser
         *   Nested (finer=inner): r_bp_query = r_outer_finer,   r_bp_candidate = ri_coarser (bore)
         *
         *
         * @param queryMid              midpoint of segment i 
         * @param queryHalfLength       half-length of segment i;
         * @param r_bp_query            contact-relevant bounding radius of the query beam
         * @param candidateFrames       Rigid3d frames of the candidate beam
         * @param r_bp_candidate        contact-relevant bounding radius of the candidate beam
         * @param broadPhaseMarginFactor  multiplier on radii (default 1.5)
         * @return                      sorted list of candidate segment indices j
         */
        static sofa::type::vector<int>
            candidateSegments(const Vec3&          queryMid,
                              Real                 queryHalfLength,
                              Real                 r_bp_query,
                              const VecRigidCoord& candidateFrames,
                              Real                 r_bp_candidate,
                              Real                 broadPhaseMarginFactor);
        
        /**
        * @brief Computes the contact-plane tangent frame {t̂₁, t̂₂} from a raw
        *        axial chord τ and the contact normal n̂.
        *
        * Algorithm (identical to SSIM evaluateContacts() inline computation):
        *   t̂₁ = normalize(τ − (τ·n̂)·n̂)      [project onto contact plane]
        *   t̂₂ = n̂ × t̂₁
        *
        * If τ ∥ n̂ (degenerate), falls back to tau2 (second chord), then global X/Y.
        *
        * @param tau1        Raw Beam-1 segment chord (unit or near-unit).
        * @param tau2        Raw Beam-2 segment chord (fallback).
        * @param nHat        Unit contact normal (raw Beam-1 → Beam-2 direction).
        * @param t1_out      [out] contact-plane axial tangent t̂₁.
        * @param t2_out      [out] contact-plane circumferential tangent t̂₂.
        */
        void computeContactFrame(const Vec3& tau1,
                                 const Vec3& tau2,
                                 const Vec3& nHat,
                                 Vec3&       t1_out,
                                 Vec3&       t2_out);
        
        /**
        * @brief Runs all parameter-consistency checks: radii positivity, inner/outer
        *        radii ordering, nested geometry validity, velocity-link sizes,
        *        defaultNormal vs contact configuration.
        *
        * Sets d_componentState to Invalid on any failure. Clears m_lastValidNormal
        * (cache may be stale after parameter changes). Called by both init() and
        * reinit() so that runtime data changes don't leave the component in a
        * stale-valid state.
        *
        * @return true if all checks pass, false otherwise.
        */
        bool validateParameters();

    private:
        static constexpr Real s_eps = Real(1e-14); ///< numerical zero

        struct EvaluationCacheKey
        {
            const sofa::core::objectmodel::BaseData* frames1 { nullptr };
            const sofa::core::objectmodel::BaseData* frames2 { nullptr };
            const sofa::core::objectmodel::BaseData* vels1 { nullptr };
            const sofa::core::objectmodel::BaseData* vels2 { nullptr };
            int frames1Counter { 0 };
            int frames2Counter { 0 };
            int vels1Counter { 0 };
            int vels2Counter { 0 };
            int parameterCounter { 0 };

            bool operator==(const EvaluationCacheKey& other) const
            {
                return frames1 == other.frames1 &&
                       frames2 == other.frames2 &&
                       vels1 == other.vels1 &&
                       vels2 == other.vels2 &&
                       frames1Counter == other.frames1Counter &&
                       frames2Counter == other.frames2Counter &&
                       vels1Counter == other.vels1Counter &&
                       vels2Counter == other.vels2Counter &&
                       parameterCounter == other.parameterCounter;
            }
        };

        /// Snapshot stored alongside each cached normal. Enables smoothness-based     
        /// validation of fallback 2a in computeContactNormal().
        /// Previously: cache stored only the Vec3 normal.
        struct CachedNormal
        {
            Vec3                   normal;   ///< cached raw Beam-1 → Beam-2 unit normal
            sofa::type::Quat<Real> qSlerp;   ///< SLERP'd Beam-1 frame orientation at contact point
        };

        /// Persistent cache of the last valid contact normal (and contact-point frame)
        /// for each section pair (i_beam1, j_beam2). Used for the zero-normal fallback.
        std::map<std::pair<int,int>, CachedNormal> m_lastValidNormal;    


        ContactEvaluation computeContacts(const VecRigidCoord& frames1,
                                          const VecRigidCoord& frames2,
                                          const VecRigidDeriv& vels1,
                                          const VecRigidDeriv& vels2);

        bool m_evalCacheValid { false };
        EvaluationCacheKey m_evalCacheKey;
        ContactEvaluation m_evalCache;
    };

} // namespace Cosserat
