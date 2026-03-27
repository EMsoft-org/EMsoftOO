# C-interop shared library for Python (and other language) bindings.
# This builds libEMsoftOO_c which exposes bind(c) wrappers around
# the EMsoftOO Fortran library.

set(C_INTERFACE_SOURCES
    c_quaternions.f90
    c_rotations.f90
    c_crystallography.f90
    c_symmetry.f90
    c_diffraction.f90
    c_so3.f90
    c_lambert.f90
)

add_library(EMsoftOO_c SHARED ${C_INTERFACE_SOURCES})

target_link_libraries(EMsoftOO_c
    EMsoftOOLib
)

target_include_directories(EMsoftOO_c
    PRIVATE ${EMsoftOOLib_BINARY_DIR}
)

set_target_properties(EMsoftOO_c PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

install(TARGETS EMsoftOO_c
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)
