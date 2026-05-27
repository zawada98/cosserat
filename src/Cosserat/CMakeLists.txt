project(Cosserat.src VERSION 21.06.99)

set(HEADER_FILES
    config.h

    mapping/BaseCosserat.h
    mapping/BaseCosserat.inl
    mapping/BeamContactMapping.h
    mapping/BeamContactMapping2.h
    mapping/DiscreteCosseratMapping.h
    mapping/DiscreteCosseratMapping.inl
    mapping/DiscretDynamicCosseratMapping.h
    mapping/DiscreteDynamicCosseratMapping.inl
    mapping/DifferenceMultiMapping.h
    mapping/DifferenceMultiMapping.inl
    mapping/RigidDistanceMapping.h
    mapping/RigidDistanceMapping.inl
    mapping/LegendrePolynomialsMapping.h
    mapping/LegendrePolynomialsMapping.inl
    
    engine/ProjectionEngine.h
    engine/ProjectionEngine.inl
    engine/ContactTriad.h



    forcefield/BeamHookeLawForceField.h
    forcefield/BeamHookeLawForceField.inl
    forcefield/BeamHookeLawForceFieldRigid.h
    forcefield/BeamHookeLawForceFieldRigid.inl
    forcefield/CosseratInternalActuation.h
    forcefield/CosseratInternalActuation.inl
    forcefield/MyUniformVelocityDampingForceField.h
    forcefield/MyUniformVelocityDampingForceField.inl

    constraint/CosseratSlidingConstraint.h
    constraint/CosseratSlidingConstraint.inl
    constraint/QPSlidingConstraint.h
    constraint/QPSlidingConstraint.inl
    constraint/CosseratNeedleSlidingConstraint.h
    constraint/CosseratNeedleSlidingConstraint.inl
    constraint/ContactPointsUnilateralConstraint.h

    intersection/SphereSweptIntersectionMethod.h
    intersection/SphereSweptIntersectionMethod2.h

    )

set(SOURCE_FILES
    initCosserat.cpp

    mapping/BaseCosserat.cpp
    mapping/BeamContactMapping.cpp
    mapping/BeamContactMapping2.cpp
    mapping/DiscreteCosseratMapping.cpp
    mapping/DiscreteDynamicCosseratMapping.cpp
    mapping/DifferenceMultiMapping.cpp
    mapping/RigidDistanceMapping.cpp
    mapping/LegendrePolynomialsMapping.cpp

    engine/ProjectionEngine.cpp

    forcefield/BeamHookeLawForceField.cpp
    forcefield/BeamHookeLawForceFieldRigid.cpp
    forcefield/CosseratInternalActuation.cpp
    forcefield/MyUniformVelocityDampingForceField.cpp

    constraint/CosseratSlidingConstraint.cpp
    constraint/QPSlidingConstraint.cpp
    constraint/CosseratNeedleSlidingConstraint.cpp
    constraint/ContactPointsUnilateralConstraint.cpp
    intersection/SphereSweptIntersectionMethod.cpp 
    )

# add_subdirectory(Binding)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>")

target_link_libraries(${PROJECT_NAME} SofaTest SofaGTestMain SofaCore
    SofaConstraint SofaBaseMechanics SofaUserInteraction)
