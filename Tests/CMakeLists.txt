cmake_minimum_required(VERSION 3.12)

set(This Cosserat_test)

project(${This} C CXX)

#find_package(Cosserat REQUIRED)
find_package(Sofa.Testing REQUIRED)
enable_testing()

set(HEADER_FILES
        Example.h
        constraint/Constraint.h
        )
set(SOURCE_FILES
#        Example.cpp
#        constraint/ExampleTest.cpp
#        constraint/CosseratUnilateralInteractionConstraintTest.cpp
 #       forcefield/BeamHookeLawForceFieldTest.cpp
	    mapping/BeamContactMappingTest.cpp
    )


add_executable(${This} ${SOURCE_FILES} ${HEADER_FILES})

target_link_libraries(${PROJECT_NAME} 
        Sofa.Testing
        Cosserat
        Sofa.Core
        Sofa.Component.Mapping.NonLinear
        C:/Users/zawada/SOFA/build/lib/Release/gtest.lib
        C:/Users/zawada/SOFA/build/lib/Release/gtest_main.lib
)

target_include_directories(${This} 
        PUBLIC 
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>"
        "C:/Users/zawada/SOFA/plugins/Cosserat/cmake-build-release-visual-studio/include/Cosserat"
        "C:/Users/zawada/SOFA/build/external_directories/fetched/googletest/googletest/include"
)


add_test(
        NAME ${This} 
        COMMAND ${This}
)

#[[add_subdirectory(constraint)]]
