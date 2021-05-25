set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/XRay")
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

endif()

if(EMsoftOO_ENABLE_HDF5_SUPPORT)

  GetHDF5LinkLibraries(EMSOFTOO)
  set(EXE_LINK_LIBRARIES ${EXE_LINK_LIBRARIES} ${EMSOFTOO_hdf5LinkLibs} )

  Add_EMsoftOO_Executable(TARGET EMLauemaster
                        SOURCES ${APP_DIR}/EMLauemaster.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMLauemaster.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )


  Add_EMsoftOO_Executable(TARGET EMLaue
                      SOURCES ${APP_DIR}/EMLaue.f90
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES}  
                      TEMPLATE ${TMPLT_DIR}/EMLaue.template
                      SOLUTION_FOLDER EMsoftPublic/XRay
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR} ${EMsoftOOHDFLib_BINARY_DIR} ${EMsoft_BINARY_DIR}/EMsoftLib
                      )

  #Add_EMsoftOO_Executable(TARGET EMLaueSlit
  #                    SOURCES ${APP_DIR}/EMLaueSlit.f90
  #                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}  
  #                    TEMPLATE ${TMPLT_DIR}/EMLaue.template
  #                    SOLUTION_FOLDER EMsoftPublic/XRay
  #                    INSTALL_PROGRAM TRUE
  #                    INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR} ${EMsoftOOHDFLib_BINARY_DIR} ${EMsoft_BINARY_DIR}/EMsoftLib
  #                   )

endif()

