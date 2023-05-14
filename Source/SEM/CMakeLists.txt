set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/SEM")
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

  Add_EMSoftOO_Executable(TARGET EMMCOpenCL
                        SOURCES ${APP_DIR}/EMMCOpenCL.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMMCOpenCL.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR} ${EMOpenCLLib_BINARY_DIR}
                        )


  # Add_EMSoftOO_Executable(TARGET EMMCfoil
  #                       SOURCES ${APP_DIR}/EMMCfoil.f90
  #                       LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                       TEMPLATE ${TMPLT_DIR}/EMMfoil.template
  #                       SOLUTION_FOLDER EMsoftOOPublic/SEM
  #                       INSTALL_PROGRAM TRUE
  #                       INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR} ${EMOpenCLLib_BINARY_DIR}
  #                       )

  # Add_EMsoftOO_Executable(TARGET EMEBSDFull
  #                       SOURCES ${APP_DIR}/EMEBSDFull.f90
  #                       LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                       TEMPLATE ${TMPLT_DIR}/EMEBSDFull.template
  #                       SOLUTION_FOLDER EMsoftOOPublic/SEM
  #                       INSTALL_PROGRAM TRUE
  #                       INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR} ${EMOpenCLLib_BINARY_DIR}
  #                       )


endif()

if(EMsoftOO_ENABLE_HDF5_SUPPORT)

  GetHDF5LinkLibraries(EMSOFTOO)
  set(EXE_LINK_LIBRARIES ${EXE_LINK_LIBRARIES} ${EMSOFTOO_hdf5LinkLibs} )

  Add_EMsoftOO_Executable(TARGET EMEBSDmaster
                        SOURCES ${APP_DIR}/EMEBSDmaster.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMEBSDmaster.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )


  Add_EMsoftOO_Executable(TARGET EMEBSD
                        SOURCES ${APP_DIR}/EMEBSD.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMEBSD.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMSoftOO_Executable(TARGET EMEBSDPCA
                        SOURCES ${APP_DIR}/EMEBSDPCA.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMEBSDPCA.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR} 
                        )

  Add_EMsoftOO_Executable(TARGET EMBSE
                        SOURCES ${APP_DIR}/EMBSE.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMBSE.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMsoftOO_Executable(TARGET EMISE
                        SOURCES ${APP_DIR}/EMISE.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMISE.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMsoftOO_Executable(TARGET EMEBSDkin
                        SOURCES ${APP_DIR}/EMEBSDkin.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMEBSDkin.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMsoftOO_Executable(TARGET EMreflectors
                        SOURCES ${APP_DIR}/EMreflectors.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMreflectors.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  # Add_EMsoftOO_Executable(TARGET EMEBSDmasterSHT
  #                       SOURCES ${APP_DIR}/EMEBSDmasterSHT.f90
  #                       LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                       TEMPLATE ${TMPLT_DIR}/EMEBSDmasterSHT.template
  #                       SOLUTION_FOLDER EMsoftOOPublic/SEM
  #                       INSTALL_PROGRAM TRUE
  #                       INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
  #                       )

  Add_EMsoftOO_Executable(TARGET EMISEmaster
                        SOURCES ${APP_DIR}/EMISEmaster.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMISEmaster.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMsoftOO_Executable(TARGET EMECP
                        SOURCES ${APP_DIR}/EMECP.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMECP.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMsoftOO_Executable(TARGET EMECPmaster
                        SOURCES ${APP_DIR}/EMECPmaster.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMECPmaster.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMsoftOO_Executable(TARGET EMKosselmaster
                        SOURCES ${APP_DIR}/EMKosselmaster.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMKosselmaster.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  # Add_EMsoftOO_Executable(TARGET EMMC
  #                       SOURCES ${APP_DIR}/EMMC.f90
  #                       LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                       TEMPLATE ${TMPLT_DIR}/EMMC.template
  #                       SOLUTION_FOLDER EMsoftOOPublic/SEM
  #                       INSTALL_PROGRAM TRUE
  #                       INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
  #                       )

   Add_EMsoftOO_Executable(TARGET EMECCI
                         SOURCES ${APP_DIR}/EMECCI.f90
                         LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                         TEMPLATE ${TMPLT_DIR}/EMECCI.template
                         SOLUTION_FOLDER EMsoftOOPublic/SEM
                         INSTALL_PROGRAM TRUE
                         INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
                         )

 #  Add_EMsoftOO_Executable(TARGET EMECCI_Bethe
 #                       SOURCES ${APP_DIR}/EMECCI_Bethe.f90
 #                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
 #                        TEMPLATE ${TMPLT_DIR}/EMECCI.template
 #                        SOLUTION_FOLDER EMsoftOOPublic/SEM
 #                        INSTALL_PROGRAM TRUE
 #                        INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
 #                        )

  # Add_EMsoftOO_Executable(TARGET EMECPSingle
  #                       SOURCES ${APP_DIR}/EMECPSingle.f90
  #                       LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                       TEMPLATE ${TMPLT_DIR}/EMECPSingle.template
  #                       SOLUTION_FOLDER EMsoftOOPublic/SEM
  #                       INSTALL_PROGRAM TRUE
  #                       INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
  #                       )

  # Add_EMsoftOO_Executable(TARGET EMKAM
  #                       SOURCES ${APP_DIR}/EMKAM.f90
  #                       LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                       TEMPLATE ${TMPLT_DIR}/EMKAM.template
  #                       SOLUTION_FOLDER EMsoftOOPublic/SEM
  #                       INSTALL_PROGRAM TRUE
  #                       INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
  #                       )

  # Add_EMsoftOO_Executable(TARGET EMAverageOrient
  #                       SOURCES ${APP_DIR}/EMAverageOrient.f90
  #                       LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                       TEMPLATE ${TMPLT_DIR}/EMAverageOrient.template
  #                       SOLUTION_FOLDER EMsoftOOPublic/SEM
  #                       INSTALL_PROGRAM TRUE
  #                       INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
  #                       )

  # Add_EMsoftOO_Executable(TARGET EMOrientationSimilarity
  #                       SOURCES ${APP_DIR}/EMOrientationSimilarity.f90
  #                       LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                       TEMPLATE ${TMPLT_DIR}/EMOrientationSimilarity.template
  #                       SOLUTION_FOLDER EMsoftOOPublic/SEM
  #                       INSTALL_PROGRAM TRUE
  #                       INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
  #                       )

  # Add_EMsoftOO_Executable(TARGET EMOrientationViz
  #                       SOURCES ${APP_DIR}/EMOrientationViz.f90
  #                       LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
  #                       TEMPLATE ${TMPLT_DIR}/EMOrientationViz.template
  #                       SOLUTION_FOLDER EMsoftOOPublic/SEM
  #                       INSTALL_PROGRAM TRUE
  #                       INCLUDE_DIRS ${EMsoftHDFLib_BINARY_DIR}
  #                       )

  Add_EMsoftOO_Executable(TARGET EMgetMP
                        SOURCES ${APP_DIR}/EMgetMP.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMgetMP.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMsoftOO_Executable(TARGET EMkinematical
                        SOURCES ${APP_DIR}/EMkinematical.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMkinematical.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMsoftOO_Executable(TARGET EMTKD
                        SOURCES ${APP_DIR}/EMTKD.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMTKD.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMsoftOO_Executable(TARGET EMTKDmaster
                        SOURCES ${APP_DIR}/EMTKDmaster.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMTKDmaster.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

  Add_EMsoftOO_Executable(TARGET EMMPoverlap
                        SOURCES ${APP_DIR}/EMMPoverlap.f90
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                        TEMPLATE ${TMPLT_DIR}/EMMPoverlap.template
                        SOLUTION_FOLDER EMsoftOOPublic/SEM
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                        )

endif()
