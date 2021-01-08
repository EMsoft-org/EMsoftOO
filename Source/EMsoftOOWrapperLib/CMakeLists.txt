
project(EMsoftOOWrapperLib)

# --------------------------------------------------------------------
# The MODALITY_DIRS variable is set one level above this one.

# --------------------------------------------------------------------
# Add the executables
foreach(MODALITY ${MODALITY_DIRS})
  if( "${EMsoftOO_ENABLE_${MODALITY}}" STREQUAL "ON" )
    # message(STATUS "${MODALITY} Wrapper Lib: ${EMsoftOO_ENABLE_${MODALITY}}")
    include( ${PROJECT_SOURCE_DIR}/${MODALITY}/SourceList.cmake)
  endif()
endforeach()


set(EMsoftOOWrapperLib_SOURCES
  ${DictionaryIndexing_Wrapper_SRCS}
  #${EEC_Wrapper_SRCS}
  #${OLIO_Wrapper_SRCS}
  ${OM_Wrapper_SRCS}
  ${SEM_Wrapper_SRCS}
  #${TEM_Wrapper_SRCS}
  #${QC_Wrapper_SRCS}
  ${Utilities_Wrapper_SRCS}
)

set(EMsoftOOWrapperLib_HEADERS
  ${DictionaryIndexing_Wrapper_HEADERS} 
  ${SEM_Wrapper_HEADERS} 
  #${TEM_Wrapper_HEADERS} 
)



add_library(EMsoftOOWrapperLib ${LIB_TYPE} ${EMsoftOOWrapperLib_SOURCES} ${EMsoftOOWrapperLib_HEADERS})
target_link_libraries(EMsoftOOWrapperLib EMsoftOOLib EMOpenCLLib)

target_include_directories(EMsoftOOWrapperLib 
      PUBLIC ${FFTW3_INCLUDE_DIR} 
      PUBLIC ${CLFORTRAN_INCLUDE_DIR}
      ${CLFortran_INSTALL}/include
      $<BUILD_INTERFACE:${EMsoftOOLib_BINARY_DIR}>
      $<BUILD_INTERFACE:${EMOpenCLLib_BINARY_DIR}>
)

set_target_properties (EMsoftOOWrapperLib PROPERTIES
  LINKER_LANGUAGE Fortran
   # Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
  FOLDER EMsoftOOWrapperPublic
)


#---------------------------------------------------------------------
# This sets up the two variables install_dir and lib_install_dir
EMsoftOO_SetupInstallDirs()

#if(BUILD_SHARED_LIBS)
  INSTALL(TARGETS EMsoftOOWrapperLib
    COMPONENT Applications
    EXPORT EMsoftOOWrapperLibTargets
    RUNTIME DESTINATION ${install_dir}
    LIBRARY DESTINATION ${lib_install_dir}
    ARCHIVE DESTINATION ${lib_install_dir}
    BUNDLE DESTINATION "."
  )
#endif()


install(FILES ${EMsoftOOWrapperLib_HEADERS}
      DESTINATION ${include_install_dir}/EMsoftOOWrapperLib
      COMPONENT Applications
)

if(APPLE AND BUILD_SHARED_LIBS)
  set(QAB_INSTALL_DEST "lib")
  set(osx_app_name "libEMsoftOOWrapperLib.dylib")
  get_property(EMsoftOOBinaryDir GLOBAL PROPERTY EMsoftOOBinaryDir)
  set(OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT "${EMsoftOO_BINARY_DIR}/EMsoftOOWrapperLib/EMsoftOOWrapperLib_CompleteTool.cmake")
  set(OPTIMIZE_BUNDLE_SHELL_SCRIPT            "${EMsoftOO_BINARY_DIR}/EMsoftOOWrapperLib/EMsoftOOWrapperLib_OptimizeTool.sh")
  set(PROJECT_INSTALL_DIR "lib")
  configure_file("${CMP_OSX_TOOLS_SOURCE_DIR}/CompleteTool.cmake.in"
                "${OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT}" @ONLY IMMEDIATE)
  configure_file("${CMP_OSX_TOOLS_SOURCE_DIR}/CompleteTool.sh.in"
                "${OPTIMIZE_BUNDLE_SHELL_SCRIPT}" @ONLY IMMEDIATE)
  #install(SCRIPT "${OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT}" COMPONENT ${QAB_COMPONENT})
endif()


# --------------------------------------------------------------------
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOWrapperLibTargetsConfigVersion.cmake"
  VERSION ${EMsoftOO_VERSION}
  COMPATIBILITY AnyNewerVersion
)
#if(BUILD_SHARED_LIBS)
  export(EXPORT EMsoftOOWrapperLibTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOWrapperLibTargets.cmake"
    NAMESPACE EMsoftOO::
  )
#endif()

configure_file(cmake/EMsoftOOWrapperConfig.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOWrapperConfig.cmake"
  @ONLY
)

set(ConfigPackageLocation ${lib_install_dir}/cmake/EMsoftOO)

#if(BUILD_SHARED_LIBS)
  install(EXPORT EMsoftOOWrapperLibTargets
    FILE
      EMsoftOOWrapperLibTargets.cmake
    NAMESPACE
      EMsoftOO::
    DESTINATION
      ${ConfigPackageLocation}
  )

#endif()

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/EMsoftOO/EMsoftOOWrapperLibTargetsConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)

