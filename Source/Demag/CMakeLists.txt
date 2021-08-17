set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/Demag")
set(TMPLT_DIR "${EMsoftOO_SOURCE_DIR}/NamelistTemplates")
set(LIB_SEARCH_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if((EMsoftOO_ENABLE_HDF5_SUPPORT) AND (EMsoftOO_ENABLE_OpenCL_SUPPORT))

  # The libraries are specifically laid out this way in order to ensure the link line
  # has the correct ordering. This seems to be more important on Windows than on macOS/Linux
  set(EXE_LINK_LIBRARIES
          ${OpenCL_LIBRARY}
          clfortran
          EMsoftOOLib
          EMOpenCLLib
  )

  # Add_EMSoftOO_Executable(TARGET EMMCOpenCL
  #                       SOURCES ${APP_DIR}/EMMCOpenCL.f90
  #                       LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                       TEMPLATE ${TMPLT_DIR}/EMMCOpenCL.template
  #                       SOLUTION_FOLDER EMsoftOOPublic/SEM
  #                       INSTALL_PROGRAM TRUE
  #                       INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR} ${EMOpenCLLib_BINARY_DIR}
  #                       )

endif()

if(EMsoftOO_ENABLE_HDF5_SUPPORT)

  GetHDF5LinkLibraries(EMSOFTOO)
  set(EXE_LINK_LIBRARIES ${EXE_LINK_LIBRARIES} ${EMSOFTOO_hdf5LinkLibs} )

  Add_EMsoftOO_Executable(TARGET EMDemag
                        SOURCES ${APP_DIR}/EMDemag.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMDemag.template
                        SOLUTION_FOLDER EMsoftOOPublic/Demag
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

endif()
