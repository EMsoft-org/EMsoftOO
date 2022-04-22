set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/DictionaryIndexing")
set(TMPLT_DIR "${EMsoftOO_SOURCE_DIR}/NamelistTemplates")
set(LIB_SEARCH_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


if((EMsoftOO_ENABLE_HDF5_SUPPORT) AND (EMsoftOO_ENABLE_OpenCL_SUPPORT))

  GetHDF5LinkLibraries(EMSOFT)
  set(EXE_LINK_LIBRARIES ${EMSOFT_hdf5LinkLibs} EMOpenCLLib)

  Add_EMsoftOO_Executable(TARGET EMDI
                      SOURCES ${APP_DIR}/EMDI.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMDI.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}  ${EMOpenCLLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMDIRAM
                      SOURCES ${APP_DIR}/EMDIRAM.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMDIRAM.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}  ${EMOpenCLLib_BINARY_DIR}
                      )

  # Add_EMsoftOO_Executable(TARGET EMECPDI
  #                     SOURCES ${APP_DIR}/EMECPDI.f90 
  #                     LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
  #                     TEMPLATE ${TMPLT_DIR}/EMECPDI.template 
  #                     SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
  #                     INSTALL_PROGRAM TRUE
  #                     INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}  ${EMOpenCLLib_BINARY_DIR}
  #                     )


  # Add_EMsoftOO_Executable(TARGET EMTKDDI
  #                     SOURCES ${APP_DIR}/EMTKDDI.f90 
  #                     LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
  #                     TEMPLATE ${TMPLT_DIR}/EMTKDDI.template 
  #                     SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
  #                     INSTALL_PROGRAM TRUE
  #                     INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}  ${EMOpenCLLib_BINARY_DIR}
  #                     )

endif()

if(EMsoftOO_ENABLE_HDF5_SUPPORT)

  GetHDF5LinkLibraries(EMSOFTOO)
  set(EXE_LINK_LIBRARIES ${EXE_LINK_LIBRARIES} ${EMSOFTOO_hdf5LinkLibs} )

  # Add_EMsoftOO_Executable(TARGET EMDPFit
  #                     SOURCES ${APP_DIR}/EMDPFit.f90 
  #                     LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                     TEMPLATE ${TMPLT_DIR}/EMDPFit.template 
  #                     SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
  #                     INSTALL_PROGRAM TRUE
  #                     INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
  #                     )

  Add_EMsoftOO_Executable(TARGET EMFitOrientation
                      SOURCES ${APP_DIR}/EMFitOrientation.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMFitOrientation.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )
 
  Add_EMsoftOO_Executable(TARGET EMgetIPF
                      SOURCES ${APP_DIR}/EMgetIPF.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMgetIPF.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMgetCTF
                      SOURCES ${APP_DIR}/EMgetCTF.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMgetCTF.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMgetANG
                      SOURCES ${APP_DIR}/EMgetANG.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMgetANG.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMgetDisOrMap
                      SOURCES ${APP_DIR}/EMgetDisOrMap.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMgetDisOrMap.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMgetEulers
                      SOURCES ${APP_DIR}/EMgetEulers.f90 
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMgetEulers.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )
 
  Add_EMsoftOO_Executable(TARGET EMDIpreview
                      SOURCES ${APP_DIR}/EMDIpreview.f90 
                      LINK_LIBRARIES jsonfortran ${EXE_LINK_LIBRARIES}  
                      TEMPLATE ${TMPLT_DIR}/EMDIpreview.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMDIchangesetting
                      SOURCES ${APP_DIR}/EMDIchangesetting.f90 
                      LINK_LIBRARIES jsonfortran ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMDIchangesetting.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMgetADP
                      SOURCES ${APP_DIR}/EMgetADP.f90 
                      LINK_LIBRARIES jsonfortran ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMgetADP.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing 
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMgetOSM
                      SOURCES ${APP_DIR}/EMgetOSM.f90 
                      LINK_LIBRARIES jsonfortran ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMgetOSM.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing 
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMdpmerge
                      SOURCES ${APP_DIR}/EMdpmerge.f90 
                      LINK_LIBRARIES jsonfortran ${EXE_LINK_LIBRARIES} 
                      TEMPLATE ${TMPLT_DIR}/EMdpmerge.template 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing 
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )

  Add_EMsoftOO_Executable(TARGET EMdpextract
                      SOURCES ${APP_DIR}/EMdpextract.f90 
                      LINK_LIBRARIES jsonfortran ${EXE_LINK_LIBRARIES} 
                      SOLUTION_FOLDER EMsoftOOPublic/DictionaryIndexing 
                      INSTALL_PROGRAM TRUE
                      INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                      )

endif()
