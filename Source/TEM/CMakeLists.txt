set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/TEM")
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


  # Add_EMsoftOO_Executable(TARGET EMdddSTEM
  #                     SOURCES ${APP_DIR}/EMdddSTEM.f90
  #                     LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                     TEMPLATE ${TMPLT_DIR}/EMmdSTEM.template 
  #                     SOLUTION_FOLDER EMsoftPublic/TEM
  #                     INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR} ${EMOpenCLLib_BINARY_DIR}
  #                     )

endif()

if(EMsoftOO_ENABLE_HDF5_SUPPORT)

  GetHDF5LinkLibraries(EMSOFTOO)
  set(EXE_LINK_LIBRARIES ${EXE_LINK_LIBRARIES} ${EMSOFTOO_hdf5LinkLibs} )

  Add_EMsoftOO_Executable(TARGET EMPEDkin 
                      SOURCES ${APP_DIR}/EMPEDkin.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMPEDkin.template
                      SOLUTION_FOLDER EMsoftOOPublic/TEM
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMLACBED
                      SOURCES ${APP_DIR}/EMLACBED.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                      TEMPLATE ${TMPLT_DIR}/EMLACBED.template 
                      SOLUTION_FOLDER EMsoftOOPublic/TEM 
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMPEDZA 
                      SOURCES ${APP_DIR}/EMPEDZA.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                      TEMPLATE ${TMPLT_DIR}/EMEBSDmaster.template 
                      SOLUTION_FOLDER EMsoftOOPublic/TEM 
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMhh4
                      SOURCES ${APP_DIR}/EMhh4.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMhh4.template
                      SOLUTION_FOLDER EMsoftOOPublic/TEM
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
                      )
 
 
endif()
