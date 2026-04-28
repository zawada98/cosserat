/******************************************************************************
 * Cosserat Plugin for SOFA Framework                                         *
 *                                                                            *
 * ContactTriad.h                                                             *
 *                                                                            *
 * Shared POD type describing one contact pair's orthonormal triad:           *
 *   n   – unit contact normal                                                *
 *   t1  – Beam-1 tangent projected onto the contact plane                   *
 *   t2  – circumferential tangent,  t2 = n × t1                             *
 *                                                                            *
 * Produced by BeamContactMapping (writes d_contactTriads per pair in         *
 * apply()), consumed by ContactPointsUnilateralConstraint (reads the triad   *
 * and uses n for the normal row, t1/t2 for the friction rows).              *
 *                                                                            *
 * ── Why a dedicated struct ───────────────────────────────────────────────── *
 *  - Keeps all three directions grouped: impossible to link "2 of 3"          *
 *    parallel Data fields in a scene file and silently get junk for the      *
 *    third axis.                                                             *
 *  - Reads like   triad.n,  triad.t1,  triad.t2 — no quaternion gymnastics,  *
 *    no hidden axis-index convention.                                        *
 *  - Zero arithmetic between producer and consumer: BCM writes three Vec3s,  *
 *    CPULC reads three Vec3s.                                                *
 *                                                                            *
 * ── DataTypeInfo support ─────────────────────────────────────────────────── *
 * SOFA's generic Data<vector<T>> serializer delegates to T's operator<< /     *
 * operator>>, so only these two are needed.  If a build ever fails on        *
 * missing DataTypeInfo, add a specialization (see SOFA doc                    *
 * "Programming with SOFA / Data in components").                             *
 ******************************************************************************/
#pragma once

#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>
#include <iostream>

namespace Cosserat
{

struct ContactTriad
{
    sofa::type::Vec3d n;   ///< unit contact normal
    sofa::type::Vec3d t1;  ///< Beam-1 tangent projected onto contact plane
    sofa::type::Vec3d t2;  ///< circumferential tangent, = n × t1

    // Equality for Data change detection / tests.
    friend bool operator==(const ContactTriad& a, const ContactTriad& b)
    {
        return a.n == b.n && a.t1 == b.t1 && a.t2 == b.t2;
    }
    friend bool operator!=(const ContactTriad& a, const ContactTriad& b)
    {
        return !(a == b);
    }

    // Stream serialization — required for Data<ContactTriad> and, via SOFA's
    // generic container serializer, for Data<vector<ContactTriad>>.
    // Format: 9 floats separated by whitespace (Vec3's operator<< emits 3).
    friend std::ostream& operator<<(std::ostream& out, const ContactTriad& c)
    {
        return out << c.n << ' ' << c.t1 << ' ' << c.t2;
    }
    friend std::istream& operator>>(std::istream& in, ContactTriad& c)
    {
        return in >> c.n >> c.t1 >> c.t2;
    }
};

using VecContactTriad = sofa::type::vector<ContactTriad>;

} // namespace Cosserat