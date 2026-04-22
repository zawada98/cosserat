# Custom FindGLEW.cmake - overrides CMake 4.2 built-in for WinDepPack
set(GLEW_INCLUDE_DIR  "C:/Users/zawada/SOFA/build/external_directories/fetched/WinDepPack/include")
set(GLEW_INCLUDE_DIRS "C:/Users/zawada/SOFA/build/external_directories/fetched/WinDepPack/include")
set(GLEW_LIBRARIES    "C:/Users/zawada/SOFA/build/external_directories/fetched/WinDepPack/lib/win64/glew32.lib")
set(GLEW_FOUND TRUE)

if(NOT TARGET GLEW::GLEW)
    add_library(GLEW::GLEW UNKNOWN IMPORTED)
    set_target_properties(GLEW::GLEW PROPERTIES
        IMPORTED_LOCATION "${GLEW_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${GLEW_INCLUDE_DIRS}"
    )
endif()