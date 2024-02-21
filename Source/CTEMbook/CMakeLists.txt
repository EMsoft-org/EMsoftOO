set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/CTEMbook")
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

  Add_EMsoftOO_Executable(TARGET EMTBSM
                      SOURCES ${APP_DIR}/EMTBSM.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      SOLUTION_FOLDER EMsoftOOPublic/CTEMbook
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMTBBFDF
                      SOURCES ${APP_DIR}/EMTBBFDF.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      SOLUTION_FOLDER EMsoftOOPublic/CTEMbook
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMSRCBED
                      SOURCES ${APP_DIR}/EMSRCBED.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMSRCBED.template 
                      SOLUTION_FOLDER EMsoftOOPublic/CTEMbook
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMBWEW
                      SOURCES ${APP_DIR}/EMBWEW.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMBWEW.template 
                      SOLUTION_FOLDER EMsoftOOPublic/CTEMbook
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
                      )

 
endif()
