set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/Utilities")
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

  Add_EMSoftOO_Executable(TARGET EMOpenCLinfo 
                        SOURCES ${APP_DIR}/EMOpenCLinfo.f90 
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                        SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR} ${EMOpenCLLib_BINARY_DIR}
                        )
endif()

if(EMsoftOO_ENABLE_HDF5_SUPPORT)

      GetHDF5LinkLibraries(EMSOFTOO)
      set(EXE_LINK_LIBRARIES ${EXE_LINK_LIBRARIES} ${EMSOFTOO_hdf5LinkLibs})
    
     #  Add_EMsoftOO_Executable(TARGET EMsoftOOSlackTest
     #                                SOURCES ${APP_DIR}/EMsoftOOSlackTest.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

      Add_EMsoftOO_Executable(TARGET EMDisorientations
                                    SOURCES ${APP_DIR}/EMDisorientations.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
                                    )

     #  Add_EMsoftOO_Executable(TARGET EMDisorientationsTwoPhase
     #                                SOURCES ${APP_DIR}/EMDisorientationsTwoPhase.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

      Add_EMsoftOO_Executable(TARGET EMConvertOrientations
                                    SOURCES ${APP_DIR}/EMConvertOrientations.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
                                    )

     #  Add_EMsoftOO_Executable(TARGET EMeqvrot
     #                                SOURCES ${APP_DIR}/EMeqvrot.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

     # Add_EMsoftOO_Executable(TARGET EMlistTC
     #                                SOURCES ${APP_DIR}/EMlistTC.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

     # Add_EMsoftOO_Executable(TARGET EMeqvPS
     #                                SOURCES ${APP_DIR}/EMeqvPS.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

     #  Add_EMsoftOO_Executable(TARGET EMqg 
     #                                SOURCES ${APP_DIR}/EMqg.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

     #  Add_EMsoftOO_Executable(TARGET EMshowxtal
     #                                SOURCES ${APP_DIR}/EMshowxtal.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

     #  Add_EMsoftOO_Executable(TARGET EMstar 
     #                                SOURCES ${APP_DIR}/EMstar.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

      Add_EMsoftOO_Executable(TARGET EMorbit 
                                    SOURCES ${APP_DIR}/EMorbit.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
                                    )

      Add_EMsoftOO_Executable(TARGET modtest
                                    SOURCES ${APP_DIR}/modtest.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
                                    )


     #  Add_EMsoftOO_Executable(TARGET EMzap 
     #                                SOURCES ${APP_DIR}/EMzap.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

     #  Add_EMsoftOO_Executable(TARGET EMZAgeom 
     #                                SOURCES ${APP_DIR}/EMZAgeom.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

     #  Add_EMsoftOO_Executable(TARGET EMorient 
     #                                SOURCES ${APP_DIR}/EMorient.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

      Add_EMsoftOO_Executable(TARGET EMlatgeom 
                                    SOURCES ${APP_DIR}/EMlatgeom.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
                                    )

     #  Add_EMsoftOO_Executable(TARGET EMstereo 
     #                                SOURCES ${APP_DIR}/EMstereo.f90 
     #                                LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
     #                                SOLUTION_FOLDER EMsoftOOPublic/Utilities 
     #                                INSTALL_PROGRAM TRUE
     #                                INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
     #                                )

      Add_EMsoftOO_Executable(TARGET EMmkxtal 
                                    SOURCES ${APP_DIR}/EMmkxtal.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                                    )

      Add_EMsoftOO_Executable(TARGET TestProgram
                                    SOURCES ${APP_DIR}/TestProgram.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                                    )

      Add_EMsoftOO_Executable(TARGET EMlistSG
                                    SOURCES ${APP_DIR}/EMlistSG.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
                                    )

      Add_EMsoftOO_Executable(TARGET EMxtalExtract 
                                    SOURCES ${APP_DIR}/EMxtalExtract.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
                                    )

      # Add_EMsoftOO_Executable(TARGET EMdrawcell 
      #                               SOURCES ${APP_DIR}/EMdrawcell.f90 
      #                               LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
      #                               SOLUTION_FOLDER EMsoftOOPublic/Utilities 
      #                               INSTALL_PROGRAM TRUE
      #                               INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
      #                               )

      # Add_EMsoftOO_Executable(TARGET EMKikuchiMap 
      #                               SOURCES ${APP_DIR}/EMKikuchiMap.f90 
      #                               LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
      #                               SOLUTION_FOLDER EMsoftOOPublic/Utilities 
      #                               INSTALL_PROGRAM TRUE
      #                               INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
      #                               )

      Add_EMsoftOO_Executable(TARGET EMfamily 
                                    SOURCES ${APP_DIR}/EMfamily.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
                                    )

      # Add_EMsoftOO_Executable(TARGET EMxtalinfo 
      #                               SOURCES ${APP_DIR}/EMxtalinfo.f90 
      #                               LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
      #                               SOLUTION_FOLDER EMsoftOOPublic/Utilities 
      #                               INSTALL_PROGRAM TRUE
      #                               INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
      #                               )

      # Add_EMsoftOO_Executable(TARGET EMHOLZ 
      #                               SOURCES ${APP_DIR}/EMHOLZ.f90 
      #                               LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
      #                               SOLUTION_FOLDER EMsoftOOPublic/Utilities 
      #                               INSTALL_PROGRAM TRUE
      #                               INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
      #                               )

      # Add_EMSoft_Executable(TARGET EMmergefiles
      #                               SOURCES ${APP_DIR}/EMmergefiles.f90 
      #                               LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
      #                               SOLUTION_FOLDER EMsoftOOPublic/Utilities 
      #                               INSTALL_PROGRAM TRUE
      #                               INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
      #                               )
      
      Add_EMSoftOO_Executable(TARGET EMsampleRFZ 
                                    SOURCES ${APP_DIR}/EMsampleRFZ.f90 
                                    LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                                    SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                                    INSTALL_PROGRAM TRUE
                                    INCLUDE_DIRS ${EMsoftOOHDFLib_BINARY_DIR}
                                    )


endif()

# Add_EMsoftOO_Executable(TARGET Qtest
#                               SOURCES ${APP_DIR}/Qtest.f90 
#                               LINK_LIBRARIES EMsoftOOLib
#                               SOLUTION_FOLDER EMsoftOOPublic/Utilities 
#                               INSTALL_PROGRAM TRUE
#                               INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
#                               )

Add_EMsoftOO_Executable(TARGET EMcuboMK
                              SOURCES ${APP_DIR}/EMcuboMK.f90 
                              LINK_LIBRARIES EMsoftOOLib
                              SOLUTION_FOLDER EMsoftOOPublic/Utilities 
                              INSTALL_PROGRAM TRUE
                              INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                              )


# Add_EMsoftOO_Executable(TARGET EMsoftOOinit 
#                               SOURCES ${APP_DIR}/EMsoftOOinit.f90 
#                               LINK_LIBRARIES EMsoftOOLib
#                               SOLUTION_FOLDER EMsoftOOPublic/Utilities 
#                               INSTALL_PROGRAM TRUE
#                               INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
#                               )

