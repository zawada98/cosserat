#include <Cosserat/intersection/SphereSweptIntersectionMethod.h>
#include <Cosserat/mapping/BeamContactMapping.h>

#include <gtest/gtest.h>

#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/accessor.h>
#include <sofa/type/Quat.h>

namespace
{
using Mapping = Cosserat::BeamContactMapping;
using Ssim = Cosserat::SphereSweptIntersectionMethod;
using RigidTypes = sofa::defaulttype::Rigid3dTypes;
using Vec3Types = sofa::defaulttype::Vec3dTypes;
using RigidObject = sofa::component::statecontainer::MechanicalObject<RigidTypes>;
using Vec3 = sofa::type::Vec3d;
using Quat = sofa::type::Quat<SReal>;

constexpr SReal tolerance = 1e-12;

void expectNearVec(const Vec3& actual, const Vec3& expected, SReal tol = tolerance)
{
    EXPECT_NEAR(actual[0], expected[0], tol);
    EXPECT_NEAR(actual[1], expected[1], tol);
    EXPECT_NEAR(actual[2], expected[2], tol);
}

Vec3 linear(const RigidTypes::Deriv& value)
{
    return value.getLinear();
}

Vec3 angular(const RigidTypes::Deriv& value)
{
    return value.getAngular();
}

RigidTypes::Coord frame(const Vec3& center)
{
    return RigidTypes::Coord(center, Quat(0, 0, 0, 1));
}

RigidTypes::Deriv deriv(const Vec3& v, const Vec3& w)
{
    RigidTypes::Deriv out;
    out.getLinear() = v;
    out.getAngular() = w;
    return out;
}

struct BeamContactMappingFixture
{
    RigidObject::SPtr beam1 = sofa::core::objectmodel::New<RigidObject>();
    RigidObject::SPtr beam2 = sofa::core::objectmodel::New<RigidObject>();
    sofa::component::statecontainer::MechanicalObject<Vec3Types>::SPtr child =
        sofa::core::objectmodel::New<sofa::component::statecontainer::MechanicalObject<Vec3Types>>();
    Ssim::SPtr ssim = sofa::core::objectmodel::New<Ssim>();
    Mapping::SPtr mapping = sofa::core::objectmodel::New<Mapping>();

    sofa::core::objectmodel::Data<Mapping::OutVecCoord> outPos;
    sofa::core::objectmodel::Data<Mapping::OutVecDeriv> outVel;
    sofa::core::objectmodel::Data<Mapping::OutVecDeriv> childForce;
    sofa::core::objectmodel::Data<Mapping::In1VecDeriv> beam1Force;
    sofa::core::objectmodel::Data<Mapping::In2VecDeriv> beam2Force;
    sofa::core::objectmodel::Data<Mapping::OutMatrixDeriv> childMatrix;
    sofa::core::objectmodel::Data<Mapping::In1MatrixDeriv> beam1Matrix;
    sofa::core::objectmodel::Data<Mapping::In2MatrixDeriv> beam2Matrix;

    BeamContactMappingFixture(bool nested = false, bool twoContacts = false)
    {
        Mapping::In1VecCoord frames1;
        Mapping::In2VecCoord frames2;
        if (nested)
        {
            frames1.push_back(frame({0.0, 0.0, 0.0}));
            frames1.push_back(frame({1.0, 0.0, 0.0}));
            frames2.push_back(frame({0.0, 0.0, 0.0}));
            frames2.push_back(frame({1.0, 0.0, 0.0}));
        }
        else if (twoContacts)
        {
            frames1.push_back(frame({0.0, 0.0, 0.0}));
            frames1.push_back(frame({1.0, 0.0, 0.0}));
            frames1.push_back(frame({2.0, 0.0, 0.0}));
            frames2.push_back(frame({0.25, 0.20, -0.50}));
            frames2.push_back(frame({0.25, 0.20, 0.50}));
            frames2.push_back(frame({1.25, 0.20, 0.50}));
        }
        else
        {
            frames1.push_back(frame({0.0, 0.0, 0.0}));
            frames1.push_back(frame({1.0, 0.0, 0.0}));
            frames2.push_back(frame({0.30, 0.20, -0.50}));
            frames2.push_back(frame({0.30, 0.20, 0.50}));
        }

        beam1->resize(frames1.size());
        beam2->resize(frames2.size());
        child->resize(0);

        {
            auto x1 = beam1->writePositions();
            auto x2 = beam2->writePositions();
            x1.wref() = frames1;
            x2.wref() = frames2;
        }

        Mapping::In1VecDeriv velocities1{
            deriv({0.10, 0.20, 0.30}, {0.00, 0.00, 0.40}),
            deriv({0.30, 0.10, 0.20}, {0.10, 0.00, 0.00}),
        };
        Mapping::In2VecDeriv velocities2{
            deriv({-0.20, 0.40, 0.10}, {0.00, 0.20, 0.00}),
            deriv({0.20, 0.30, -0.10}, {0.00, 0.00, -0.30}),
        };
        if (twoContacts)
        {
            velocities1.push_back(deriv({0.05, -0.10, 0.40}, {0.20, -0.10, 0.10}));
            velocities2.push_back(deriv({-0.10, 0.15, 0.25}, {-0.20, 0.10, 0.30}));
        }

        {
            auto v1 = beam1->writeVelocities();
            auto v2 = beam2->writeVelocities();
            v1.wref() = velocities1;
            v2.wref() = velocities2;
        }

        ssim->d_beam1Frames.setParent(beam1->findData("position"));
        ssim->d_beam2Frames.setParent(beam2->findData("position"));
        ssim->d_beam1Velocities.setParent(beam1->findData("velocity"));
        ssim->d_beam2Velocities.setParent(beam2->findData("velocity"));
        ssim->d_radius1.setValue(nested ? 0.20 : 0.10);
        ssim->d_radius2.setValue(0.10);
        ssim->d_innerRadius1.setValue(nested ? 0.15 : 0.0);
        ssim->d_innerRadius2.setValue(0.0);
        ssim->d_defaultNormal.setValue({0.0, 1.0, 0.0});

        auto configuration = ssim->d_contactConfiguration.getValue();
        configuration.setSelectedItem(nested ? "nested" : "external");
        ssim->d_contactConfiguration.setValue(configuration);
        ssim->init();

        mapping->l_ssim.set(ssim.get());
        mapping->addInputModel1(beam1.get());
        mapping->addInputModel2(beam2.get());
        mapping->addOutputModel(child.get());
        mapping->init();

        beam1Force.setValue(Mapping::In1VecDeriv(frames1.size()));
        beam2Force.setValue(Mapping::In2VecDeriv(frames2.size()));
    }

    void setMode(const std::string& mode)
    {
        mapping->d_mappingMode.setValue(mode);
        mapping->reinit();
    }

    void apply()
    {
        const auto* x1 = beam1->read(sofa::core::vec_id::read_access::position);
        const auto* x2 = beam2->read(sofa::core::vec_id::read_access::position);

        const sofa::type::vector<sofa::core::objectmodel::Data<Mapping::OutVecCoord>*> out{&outPos};
        const sofa::type::vector<const sofa::core::objectmodel::Data<Mapping::In1VecCoord>*> in1{x1};
        const sofa::type::vector<const sofa::core::objectmodel::Data<Mapping::In2VecCoord>*> in2{x2};

        mapping->apply(sofa::core::mechanicalparams::defaultInstance(), out, in1, in2);
    }

    void applyJ()
    {
        const auto* v1 = beam1->read(sofa::core::vec_id::read_access::velocity);
        const auto* v2 = beam2->read(sofa::core::vec_id::read_access::velocity);

        const sofa::type::vector<sofa::core::objectmodel::Data<Mapping::OutVecDeriv>*> out{&outVel};
        const sofa::type::vector<const sofa::core::objectmodel::Data<Mapping::In1VecDeriv>*> in1{v1};
        const sofa::type::vector<const sofa::core::objectmodel::Data<Mapping::In2VecDeriv>*> in2{v2};

        mapping->applyJ(sofa::core::mechanicalparams::defaultInstance(), out, in1, in2);
    }

    void applyJT()
    {
        const sofa::type::vector<sofa::core::objectmodel::Data<Mapping::In1VecDeriv>*> out1{&beam1Force};
        const sofa::type::vector<sofa::core::objectmodel::Data<Mapping::In2VecDeriv>*> out2{&beam2Force};
        const sofa::type::vector<const sofa::core::objectmodel::Data<Mapping::OutVecDeriv>*> in{&childForce};

        mapping->applyJT(sofa::core::mechanicalparams::defaultInstance(), out1, out2, in);
    }

    void applyJTMatrix()
    {
        const sofa::type::vector<sofa::core::objectmodel::Data<Mapping::In1MatrixDeriv>*> out1{&beam1Matrix};
        const sofa::type::vector<sofa::core::objectmodel::Data<Mapping::In2MatrixDeriv>*> out2{&beam2Matrix};
        const sofa::type::vector<const sofa::core::objectmodel::Data<Mapping::OutMatrixDeriv>*> in{&childMatrix};

        mapping->applyJT(sofa::core::constraintparams::defaultInstance(), out1, out2, in);
    }

    Ssim::ContactEvaluation evaluateContacts()
    {
        return ssim->evaluateContacts(*beam1->read(sofa::core::vec_id::read_access::position),
                                            *beam2->read(sofa::core::vec_id::read_access::position),
                                            *beam1->read(sofa::core::vec_id::read_access::velocity),
                                            *beam2->read(sofa::core::vec_id::read_access::velocity));
    }
};

std::size_t contactCount(BeamContactMappingFixture& fixture)
{
    return fixture.evaluateContacts().contactNormals.size();
}

Vec3 contactPointVelocity(BeamContactMappingFixture& fixture,
                          std::size_t k,
                          bool beam1)
{
    const auto eval = fixture.evaluateContacts();
    const auto sections = eval.contactSectionIds[k];
    const auto params = eval.curvilinearParams[k];
    const auto& velocities = beam1 ? fixture.beam1->readVelocities().ref()
                                   : fixture.beam2->readVelocities().ref();
    const int i = beam1 ? sections[0] : sections[1];
    const SReal alpha = beam1 ? params[0] : params[1];
    const Vec3 contactPoint = beam1 ? eval.surfacePoints1[k] : eval.surfacePoints2[k];
    const Vec3 centerlinePoint = beam1 ? eval.centerlinePoints1[k] : eval.centerlinePoints2[k];
    const Vec3 arm = contactPoint - centerlinePoint;

    return (velocities[i].getLinear() + sofa::type::cross(velocities[i].getAngular(), arm)) * (1.0 - alpha)
         + (velocities[i + 1].getLinear() + sofa::type::cross(velocities[i + 1].getAngular(), arm)) * alpha;
}

void addPointForce(Mapping::In1VecDeriv& forces,
                   int frameIndex,
                   SReal weight,
                   const Vec3& arm,
                   const Vec3& force,
                   SReal sign)
{
    forces[frameIndex].getLinear() += force * (sign * weight);
    forces[frameIndex].getAngular() += sofa::type::cross(arm, force) * (sign * weight);
}

Mapping::In1VecDeriv expectedBackProjection(BeamContactMappingFixture& fixture,
                                             std::size_t k,
                                             const Vec3& force,
                                             bool beam1,
                                             SReal sign,
                                             std::size_t size)
{
    const auto eval = fixture.evaluateContacts();
    Mapping::In1VecDeriv out(size);
    const auto sections = eval.contactSectionIds[k];
    const auto params = eval.curvilinearParams[k];
    const int i = beam1 ? sections[0] : sections[1];
    const SReal alpha = beam1 ? params[0] : params[1];
    const Vec3 contactPoint = beam1 ? eval.surfacePoints1[k] : eval.surfacePoints2[k];
    const Vec3 centerlinePoint = beam1 ? eval.centerlinePoints1[k] : eval.centerlinePoints2[k];
    const Vec3 arm = contactPoint - centerlinePoint;

    addPointForce(out, i, 1.0 - alpha, arm, force, sign);
    addPointForce(out, i + 1, alpha, arm, force, sign);
    return out;
}

template<class Matrix>
typename Matrix::Data getMatrixCol(const Matrix& matrix, unsigned int row, unsigned int col)
{
    for (auto rowIt = matrix.begin(); rowIt != matrix.end(); ++rowIt)
    {
        if (rowIt.index() != row)
            continue;
        for (auto colIt = rowIt.begin(); colIt != rowIt.end(); ++colIt)
        {
            if (colIt.index() == col)
                return colIt.val();
        }
    }
    return typename Matrix::Data{};
}
} // namespace

// Verifies the real BCM apply() implementation in contact-points mode.
// The output must use one interleaved Vec3 MechanicalObject:
//   out[2k]   = beam-1 surface contact point from SSIM
//   out[2k+1] = beam-2 surface contact point from SSIM
TEST(BeamContactMapping, ApplyContactPointsModeCallsRealMappingAndUsesEvenOddLayout)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    fixture.setMode("contactPoints");
    fixture.apply();

    const auto& out = fixture.outPos.getValue();
    ASSERT_EQ(out.size(), 2u);
    expectNearVec(out[0], fixture.evaluateContacts().surfacePoints1[0]);
    expectNearVec(out[1], fixture.evaluateContacts().surfacePoints2[0]);
}

// Verifies the real BCM apply() implementation in gap mode.
// BCM must copy SSIM's local gap vector to its child output and publish the
// same per-contact distance/triad Data consumed by downstream contact code.
TEST(BeamContactMapping, ApplyGapModeCallsRealMappingAndPublishesSsimGapAndTriad)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    fixture.setMode("gap");
    fixture.apply();

    const auto& out = fixture.outPos.getValue();
    ASSERT_EQ(out.size(), 1u);
    expectNearVec(out[0], fixture.evaluateContacts().distances[0]);

    const auto& distances = fixture.mapping->d_distances.getValue();
    const auto& triads = fixture.mapping->d_contactTriads.getValue();
    ASSERT_EQ(distances.size(), 1u);
    ASSERT_EQ(triads.size(), 1u);
    expectNearVec(distances[0], fixture.evaluateContacts().distances[0]);
    expectNearVec(triads[0].n, fixture.evaluateContacts().contactNormals[0]);
    expectNearVec(triads[0].t1, fixture.evaluateContacts().contactTangents1[0]);
    expectNearVec(triads[0].t2, fixture.evaluateContacts().contactTangents2[0]);
}

// Verifies that the fixture contact is intentionally not an endpoint-only case.
// This protects the tests from passing while BCM ignores one interpolation
// frame, since both alpha and beta must contribute to the Jacobian cache.
TEST(BeamContactMapping, SsimFixtureUsesInteriorInterpolationWeights)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    const auto params = fixture.evaluateContacts().curvilinearParams[0];
    EXPECT_GT(params[0], 0.0);
    EXPECT_LT(params[0], 1.0);
    EXPECT_GT(params[1], 0.0);
    EXPECT_LT(params[1], 1.0);
}

// Verifies the real BCM applyJ() implementation in contact-points mode.
// The even child velocity must be the beam-1 surface-point velocity and the odd
// child velocity must be the beam-2 surface-point velocity, including angular
// velocity contributions and interpolation across both bounding frames.
TEST(BeamContactMapping, ApplyJContactPointsModeMapsBothSurfacePointVelocities)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    fixture.setMode("contactPoints");
    fixture.applyJ();

    const auto& out = fixture.outVel.getValue();
    ASSERT_EQ(out.size(), 2u);
    expectNearVec(out[0], contactPointVelocity(fixture, 0, true));
    expectNearVec(out[1], contactPointVelocity(fixture, 0, false));
}

// Verifies the real BCM applyJ() implementation in gap mode.
// Parent rigid-frame velocities must be mapped to the relative contact-point
// velocity and projected onto SSIM's {normal, tangent1, tangent2} basis.
TEST(BeamContactMapping, ApplyJGapModeProjectsRelativeVelocityIntoContactBasis)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    fixture.setMode("gap");
    fixture.applyJ();

    const auto& out = fixture.outVel.getValue();
    ASSERT_EQ(out.size(), 1u);

    const Vec3 normal = fixture.evaluateContacts().contactNormals[0];
    const Vec3 t1 = fixture.evaluateContacts().contactTangents1[0];
    const Vec3 t2 = fixture.evaluateContacts().contactTangents2[0];
    const SReal sign = fixture.ssim->gapSignForPublishedNormal();

    const Vec3 vcA = contactPointVelocity(fixture, 0, true);
    const Vec3 vcB = contactPointVelocity(fixture, 0, false);
    const Vec3 dv = vcB - vcA;

    expectNearVec(out[0], {sign * (dv * normal), dv * t1, dv * t2});
}

// applyJ may evaluate SSIM for velocity propagation, but it must not publish
// contact data or replace the frozen cache used by applyJT.  That cache belongs
// to apply(), matching SOFA's corrected-position propagation order.
TEST(BeamContactMapping, ApplyJDoesNotPublishOrOverwriteFrozenConstraintCache)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    fixture.setMode("contactPoints");
    fixture.apply();
    const auto frozenTriads = fixture.mapping->d_contactTriads.getValue();
    const auto frozenDistances = fixture.mapping->d_distances.getValue();

    fixture.childForce.setValue(Mapping::OutVecDeriv{
        Vec3(2.0, 0.0, 0.0),
        Vec3(0.0, 3.0, 0.0),
    });
    fixture.applyJT();
    const auto frozenForce1 = fixture.beam1Force.getValue();
    const auto frozenForce2 = fixture.beam2Force.getValue();

    {
        auto x2 = fixture.beam2->writePositions();
        x2.wref()[0] = frame({0.30, 1.20, -0.50});
        x2.wref()[1] = frame({0.30, 1.20,  0.50});
    }

    fixture.applyJ();
    EXPECT_EQ(fixture.mapping->d_contactTriads.getValue().size(), frozenTriads.size());
    EXPECT_EQ(fixture.mapping->d_distances.getValue().size(), frozenDistances.size());
    for (std::size_t k = 0; k < frozenTriads.size(); ++k)
    {
        expectNearVec(fixture.mapping->d_contactTriads.getValue()[k].n, frozenTriads[k].n);
        expectNearVec(fixture.mapping->d_contactTriads.getValue()[k].t1, frozenTriads[k].t1);
        expectNearVec(fixture.mapping->d_contactTriads.getValue()[k].t2, frozenTriads[k].t2);
        expectNearVec(fixture.mapping->d_distances.getValue()[k], frozenDistances[k]);
    }

    fixture.beam1Force.setValue(Mapping::In1VecDeriv(frozenForce1.size()));
    fixture.beam2Force.setValue(Mapping::In2VecDeriv(frozenForce2.size()));
    fixture.applyJT();

    const auto& force1 = fixture.beam1Force.getValue();
    const auto& force2 = fixture.beam2Force.getValue();
    ASSERT_EQ(force1.size(), frozenForce1.size());
    ASSERT_EQ(force2.size(), frozenForce2.size());
    for (std::size_t i = 0; i < force1.size(); ++i)
    {
        expectNearVec(linear(force1[i]), linear(frozenForce1[i]));
        expectNearVec(angular(force1[i]), angular(frozenForce1[i]));
    }
    for (std::size_t i = 0; i < force2.size(); ++i)
    {
        expectNearVec(linear(force2[i]), linear(frozenForce2[i]));
        expectNearVec(angular(force2[i]), angular(frozenForce2[i]));
    }
}

// Verifies the real BCM applyJT(VecDeriv) implementation in contact-points mode.
// Even child force slots must back-project only to beam 1, and odd child force
// slots must back-project only to beam 2.
TEST(BeamContactMapping, ApplyJTContactPointsModeBackProjectsEvenAndOddForcesToSeparateBeams)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    fixture.setMode("contactPoints");
    fixture.childForce.setValue(Mapping::OutVecDeriv{
        Vec3(2.0, 0.0, 0.0),
        Vec3(0.0, 3.0, 0.0),
    });
    fixture.applyJT();

    const auto& f1 = fixture.beam1Force.getValue();
    const auto& f2 = fixture.beam2Force.getValue();
    ASSERT_EQ(f1.size(), 2u);
    ASSERT_EQ(f2.size(), 2u);

    const auto expected1 = expectedBackProjection(fixture, 0, {2.0, 0.0, 0.0}, true, 1.0, f1.size());
    const auto expected2 = expectedBackProjection(fixture, 0, {0.0, 3.0, 0.0}, false, 1.0, f2.size());
    for (std::size_t i = 0; i < f1.size(); ++i)
    {
        expectNearVec(linear(f1[i]), linear(expected1[i]));
        expectNearVec(angular(f1[i]), angular(expected1[i]));
    }
    for (std::size_t i = 0; i < f2.size(); ++i)
    {
        expectNearVec(linear(f2[i]), linear(expected2[i]));
        expectNearVec(angular(f2[i]), angular(expected2[i]));
    }
}

// Verifies the real BCM applyJT(VecDeriv) implementation in gap mode.
// Local gap-space force components must be converted to world space with
// gapSign and contact tangents, then applied with opposite signs to both beams.
TEST(BeamContactMapping, ApplyJTGapModeConvertsLocalForceAndUsesOppositeBeamSigns)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    fixture.setMode("gap");
    fixture.childForce.setValue(Mapping::OutVecDeriv{Vec3(1.5, 0.25, -0.50)});
    fixture.applyJT();

    const Vec3 normal = fixture.evaluateContacts().contactNormals[0];
    const Vec3 t1 = fixture.evaluateContacts().contactTangents1[0];
    const Vec3 t2 = fixture.evaluateContacts().contactTangents2[0];
    const SReal sign = fixture.ssim->gapSignForPublishedNormal();
    const Vec3 worldForce = sign * normal * 1.5 + t1 * 0.25 + t2 * -0.50;

    const auto& f1 = fixture.beam1Force.getValue();
    const auto& f2 = fixture.beam2Force.getValue();

    const auto expected1 = expectedBackProjection(fixture, 0, worldForce, true, -1.0, f1.size());
    const auto expected2 = expectedBackProjection(fixture, 0, worldForce, false, 1.0, f2.size());
    for (std::size_t i = 0; i < f1.size(); ++i)
    {
        expectNearVec(linear(f1[i]), linear(expected1[i]));
        expectNearVec(angular(f1[i]), angular(expected1[i]));
    }
    for (std::size_t i = 0; i < f2.size(); ++i)
    {
        expectNearVec(linear(f2[i]), linear(expected2[i]));
        expectNearVec(angular(f2[i]), angular(expected2[i]));
    }
}

// Verifies the real BCM applyJT(MatrixDeriv) implementation in gap mode.
// A matrix row written in contact-local coordinates must produce the same
// parent Jacobian blocks as applyJT(VecDeriv) would produce parent forces.
TEST(BeamContactMapping, MatrixApplyJTGapModeMatchesVectorApplyJT)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    fixture.setMode("gap");
    const Vec3 localDirection(1.5, 0.25, -0.50);
    {
        auto* matrix = fixture.childMatrix.beginEdit();
        auto row = matrix->writeLine(7);
        row.addCol(0, localDirection);
        fixture.childMatrix.endEdit();
    }

    fixture.childForce.setValue(Mapping::OutVecDeriv{localDirection});
    fixture.applyJT();
    fixture.applyJTMatrix();

    const auto& vector1 = fixture.beam1Force.getValue();
    const auto& vector2 = fixture.beam2Force.getValue();
    const auto& matrix1 = fixture.beam1Matrix.getValue();
    const auto& matrix2 = fixture.beam2Matrix.getValue();
    for (std::size_t i = 0; i < vector1.size(); ++i)
        expectNearVec(getMatrixCol<Mapping::In1MatrixDeriv>(matrix1, 7, static_cast<unsigned int>(i)).getLinear(),
                      linear(vector1[i]));
    for (std::size_t i = 0; i < vector2.size(); ++i)
        expectNearVec(getMatrixCol<Mapping::In2MatrixDeriv>(matrix2, 7, static_cast<unsigned int>(i)).getLinear(),
                      linear(vector2[i]));
}

// Verifies the real BCM applyJT(MatrixDeriv) implementation in contact-points mode.
// Matrix columns must follow the same even/odd routing as vector forces:
// even columns affect only beam 1, odd columns affect only beam 2.
TEST(BeamContactMapping, MatrixApplyJTContactPointsModeUsesEvenOddRouting)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    fixture.setMode("contactPoints");
    {
        auto* matrix = fixture.childMatrix.beginEdit();
        auto row = matrix->writeLine(11);
        row.addCol(0, Vec3(2.0, 0.0, 0.0));
        row.addCol(1, Vec3(0.0, 3.0, 0.0));
        fixture.childMatrix.endEdit();
    }
    fixture.applyJTMatrix();

    const auto expected1 = expectedBackProjection(fixture, 0, {2.0, 0.0, 0.0}, true, 1.0, 2);
    const auto expected2 = expectedBackProjection(fixture, 0, {0.0, 3.0, 0.0}, false, 1.0, 2);
    const auto& matrix1 = fixture.beam1Matrix.getValue();
    const auto& matrix2 = fixture.beam2Matrix.getValue();
    for (unsigned int i = 0; i < 2; ++i)
    {
        expectNearVec(getMatrixCol<Mapping::In1MatrixDeriv>(matrix1, 11, i).getLinear(), linear(expected1[i]));
        expectNearVec(getMatrixCol<Mapping::In2MatrixDeriv>(matrix2, 11, i).getLinear(), linear(expected2[i]));
    }
}

// Verifies BCM handles more than one SSIM contact in a single apply/applyJ pass.
// The output must grow to the per-mode contact count and each slot must match
// the corresponding SSIM contact.
TEST(BeamContactMapping, MultipleContactsAreMappedIndependently)
{
    BeamContactMappingFixture fixture(false, true);
    ASSERT_GE(contactCount(fixture), 2u);

    fixture.setMode("contactPoints");
    fixture.apply();
    fixture.applyJ();

    const auto& out = fixture.outPos.getValue();
    const auto& vel = fixture.outVel.getValue();
    ASSERT_EQ(out.size(), 2 * contactCount(fixture));
    ASSERT_EQ(vel.size(), 2 * contactCount(fixture));
    for (std::size_t k = 0; k < contactCount(fixture); ++k)
    {
        const auto eval = fixture.evaluateContacts();
        expectNearVec(out[2 * k], eval.surfacePoints1[k]);
        expectNearVec(out[2 * k + 1], eval.surfacePoints2[k]);
        expectNearVec(vel[2 * k], contactPointVelocity(fixture, k, true));
        expectNearVec(vel[2 * k + 1], contactPointVelocity(fixture, k, false));
    }
}

// Verifies nested contact sign handling through the real SSIM/BCM link.
// For Beam 1 as the outer tube, SSIM publishes gapSign=-1 and BCM must use it
// consistently in gap-mode applyJ and applyJT.
TEST(BeamContactMapping, NestedGapModeUsesPublishedNegativeGapSign)
{
    BeamContactMappingFixture fixture(true);
    ASSERT_EQ(contactCount(fixture), 1u);
    ASSERT_NEAR(fixture.ssim->gapSignForPublishedNormal(), -1.0, tolerance);

    fixture.setMode("gap");
    fixture.applyJ();

    const Vec3 vcA = contactPointVelocity(fixture, 0, true);
    const Vec3 vcB = contactPointVelocity(fixture, 0, false);
    const Vec3 dv = vcB - vcA;
    const Vec3 expected(fixture.ssim->gapSignForPublishedNormal() * (dv * fixture.evaluateContacts().contactNormals[0]),
                        dv * fixture.evaluateContacts().contactTangents1[0],
                        dv * fixture.evaluateContacts().contactTangents2[0]);
    ASSERT_EQ(fixture.outVel.getValue().size(), 1u);
    expectNearVec(fixture.outVel.getValue()[0], expected);
}

// Verifies BCM's monotonic resize behavior in apply().
// Existing child slots beyond the active contact count must be preserved so
// stale constraint indices cannot become out-of-bounds during solver write-back.
TEST(BeamContactMapping, ApplyDoesNotShrinkExistingOutputSlots)
{
    BeamContactMappingFixture fixture;
    ASSERT_EQ(contactCount(fixture), 1u);

    fixture.setMode("contactPoints");
    fixture.outPos.setValue(Mapping::OutVecCoord{
        Vec3(9.0, 9.0, 9.0),
        Vec3(8.0, 8.0, 8.0),
        Vec3(7.0, 7.0, 7.0),
        Vec3(6.0, 6.0, 6.0),
    });

    fixture.apply();

    const auto& out = fixture.outPos.getValue();
    ASSERT_EQ(out.size(), 4u);
    expectNearVec(out[0], fixture.evaluateContacts().surfacePoints1[0]);
    expectNearVec(out[1], fixture.evaluateContacts().surfacePoints2[0]);
    expectNearVec(out[2], {7.0, 7.0, 7.0});
    expectNearVec(out[3], {6.0, 6.0, 6.0});
}
