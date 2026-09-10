include(${PROJECT_SOURCE_DIR}/Source/EMsoftOO_Functions.cmake)

project(EMOpenCLLib)

set(install_dir "bin")
set(lib_install_dir "lib")

if(APPLE)
  get_property(EMsoftOO_PACKAGE_DEST_PREFIX GLOBAL PROPERTY EMsoftOO_PACKAGE_DEST_PREFIX)
  set(install_dir "${EMsoftOO_PACKAGE_DEST_PREFIX}/bin")
  set(lib_install_dir "${EMsoftOO_PACKAGE_DEST_PREFIX}/lib")
elseif(WIN32)
  set(install_dir ".")
  set(lib_install_dir ".")
endif()


# GPU backend source: the OpenCL mod_GPUsupport.f90, OR (when Metal is enabled)
# the drop-in mod_GPUsupport_metal.f90 + the metal-cpp C shim.  Both provide the
# identical mod_GPUsupport / GPU_T interface, so the program modules below are
# unchanged.  See MetalMigrationPlan.md (Phase 1).
if(EMsoftOO_ENABLE_Metal_SUPPORT)
  set(EMOpenCLLib_GPU_BACKEND_SRCS
    ${EMOpenCLLib_SOURCE_DIR}/clfortran_metal_stub.f90   # provides 'use clfortran' constants
    ${EMOpenCLLib_SOURCE_DIR}/mod_GPUsupport_metal.f90
    ${EMOpenCLLib_SOURCE_DIR}/metal/emtl_shim.cpp
  )
else()
  set(EMOpenCLLib_GPU_BACKEND_SRCS
    ${EMOpenCLLib_SOURCE_DIR}/mod_GPUsupport.f90
  )
endif()

set(EMOpenCLLib_SRCS
  ${EMOpenCLLib_GPU_BACKEND_SRCS}
  # ${EMOpenCLLib_SOURCE_DIR}/Indexingmod.f90

  ${EMOpenCLLib_SOURCE_DIR}/program_mods/mod_MCOpenCL.f90
  ${EMOpenCLLib_SOURCE_DIR}/program_mods/mod_EBSDFull.f90
  ${EMOpenCLLib_SOURCE_DIR}/program_mods/mod_SEMCLwrappers.f90
  ${EMOpenCLLib_SOURCE_DIR}/program_mods/mod_DI.f90
  # ${EMOpenCLLib_SOURCE_DIR}/program_mods/mod_DIPCA.f90
  ${EMOpenCLLib_SOURCE_DIR}/program_mods/mod_HROSM.f90

# 	${EMOpenCLLib_SOURCE_DIR}/EMdymodHDF.f90
)

GetHDF5LinkLibraries(EMSOFTOO)

add_library(EMOpenCLLib ${LIB_TYPE} ${EMOpenCLLib_SRCS})

if(EMsoftOO_ENABLE_Metal_SUPPORT)
  # ---- Metal backend: metal-cpp (vendored) + Apple frameworks ----
  find_library(METAL_FRAMEWORK      Metal      REQUIRED)
  find_library(FOUNDATION_FRAMEWORK Foundation REQUIRED)
  find_library(QUARTZCORE_FRAMEWORK QuartzCore REQUIRED)

  set_source_files_properties(${EMOpenCLLib_SOURCE_DIR}/metal/emtl_shim.cpp
                              PROPERTIES COMPILE_FLAGS "-std=c++17")

  target_include_directories(EMOpenCLLib
                            PUBLIC
                              $<BUILD_INTERFACE:${EMsoftOOLib_BINARY_DIR}>
                            PRIVATE
                              "${EMsoftOO_SOURCE_DIR}/Source"
                              "${EMsoftOO_SOURCE_DIR}/ExternalProjects/metal-cpp"
                              "${EMOpenCLLib_SOURCE_DIR}/metal"
                      )

  target_link_libraries(EMOpenCLLib
                            ${METAL_FRAMEWORK}
                            ${FOUNDATION_FRAMEWORK}
                            ${QUARTZCORE_FRAMEWORK}
                            EMsoftOOLib
  )
else()
  # ---- OpenCL backend (clfortran + OpenCL runtime) ----
  target_include_directories(EMOpenCLLib
                            PUBLIC
                              ${CLFortran_INSTALL}/include
                              $<BUILD_INTERFACE:${EMsoftOOLib_BINARY_DIR}>
                            PRIVATE
                              "${EMsoftOO_SOURCE_DIR}/Source"
                      )

  target_link_libraries(EMOpenCLLib
                            ${OpenCL_LIBRARY}
                            clfortran
                            EMsoftOOLib
  )
endif()


if(EMsoftOO_ENABLE_Metal_SUPPORT)
  # libc++ is needed at final link because the shim uses the C++ standard library
  target_link_libraries(EMOpenCLLib c++)

  # Compile the .metal kernels to .metallib at build time (xcrun metal/metallib)
  # and place them next to the .cl files in the runtime opencl folder, so that
  # mod_GPUsupport_metal's read_source_file resolves them via OpenCLpathname.
  set(EMsoftOO_METAL_KERNELS
    EMMC              # Phase 1 (MC: full/bse1 modes)
    DictIndx          # Phase 2 (InnerProd)
    EMMCfoil          # Phase 3 (MC: foil mode)
    EMMCxyz           # Phase 3 (MCxyz: Ivol mode)
    # MBmoduleOpenCL is dead code (no compiled module uses it; master patterns
    # use the CPU ZGEEV CalcLgh in mod_gvectors) -- intentionally not built.
  )
  if(NOT EXISTS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl")
    file(MAKE_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl")
  endif()
  # The .metallib must sit beside the .cl files wherever OpenCLpathname resolves.
  # We emit it into the build's Bin/opencl/ (packaged/installed runs) AND copy it
  # into the source-tree opencl/ folder (the typical EMdevelop runtime location,
  # where the .cl files are version-controlled). opencl/*.metallib is gitignored.
  set(EMsoftOO_METALLIB_OUTPUTS "")
  foreach(k ${EMsoftOO_METAL_KERNELS})
    add_custom_command(
      OUTPUT  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl/${k}.metallib
      COMMAND xcrun -sdk macosx metal    -c ${EMOpenCLLib_SOURCE_DIR}/metal/${k}.metal
                                          -o ${CMAKE_CURRENT_BINARY_DIR}/${k}.air
      COMMAND xcrun -sdk macosx metallib    ${CMAKE_CURRENT_BINARY_DIR}/${k}.air
                                          -o ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl/${k}.metallib
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl/${k}.metallib
                                       ${EMsoftOO_SOURCE_DIR}/opencl/${k}.metallib
      DEPENDS ${EMOpenCLLib_SOURCE_DIR}/metal/${k}.metal
      COMMENT "Building Metal library ${k}.metallib"
      VERBATIM
    )
    list(APPEND EMsoftOO_METALLIB_OUTPUTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl/${k}.metallib)
  endforeach()
  add_custom_target(EMsoftOO_metallibs ALL DEPENDS ${EMsoftOO_METALLIB_OUTPUTS})
  add_dependencies(EMOpenCLLib EMsoftOO_metallibs)
endif()

set_target_properties (EMOpenCLLib PROPERTIES
  LINKER_LANGUAGE Fortran
  FOLDER EMsoftOOPublic
  BUILD_RPATH "${EMsoftOO_OpenMP_LIB_DIR}"
)

EMsoftOO_SetupInstallDirs()

#if(BUILD_SHARED_LIBS)
    INSTALL(TARGETS EMOpenCLLib
      EXPORT EMOpenCLLibTargets
      COMPONENT Applications
      RUNTIME DESTINATION ${install_dir}
      LIBRARY DESTINATION ${lib_install_dir}
      ARCHIVE DESTINATION ${lib_install_dir}
      BUNDLE DESTINATION "."
  )
#endif()

if(APPLE AND BUILD_SHARED_LIBS)
  set(QAB_INSTALL_DEST "lib")
  set(osx_app_name "libEMOpenCLLib.dylib")
  get_property(EMsoftOOBinaryDir GLOBAL PROPERTY EMsoftOOBinaryDir)
  set(OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT "${EMsoftOO_BINARY_DIR}/EMsoftOOLib/EMsoftOOLib_CompleteTool.cmake")
  set(OPTIMIZE_BUNDLE_SHELL_SCRIPT            "${EMsoftOO_BINARY_DIR}/EMsoftOOLib/EMsoftOOLib_OptimizeTool.sh")
  set(PROJECT_INSTALL_DIR "lib")
  configure_file("${CMP_OSX_TOOLS_SOURCE_DIR}/CompleteTool.cmake.in"
                "${OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT}" @ONLY IMMEDIATE)
  configure_file("${CMP_OSX_TOOLS_SOURCE_DIR}/CompleteTool.sh.in"
                "${OPTIMIZE_BUNDLE_SHELL_SCRIPT}" @ONLY IMMEDIATE)
  #install(SCRIPT "${OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT}" COMPONENT ${QAB_COMPONENT})
endif()

# --------------------------------------------------------------------
#
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMOpenCLLibTargetsConfigVersion.cmake"
  VERSION ${EMsoftOO_VERSION}
  COMPATIBILITY AnyNewerVersion
)
#if(BUILD_SHARED_LIBS)
  export(EXPORT EMOpenCLLibTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMOpenCLLibTargets.cmake"
    NAMESPACE EMsoftOO::
  )
#endif()

set(ConfigPackageLocation ${lib_install_dir}/cmake/EMsoftOO)

#if(BUILD_SHARED_LIBS)
  install(EXPORT EMOpenCLLibTargets
    FILE
      EMOpenCLLibTargets.cmake
    NAMESPACE
      EMsoftOO::
    DESTINATION
      ${ConfigPackageLocation}
  )
#endif()
