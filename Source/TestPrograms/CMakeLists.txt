set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/TestPrograms")
set(TMPLT_DIR "${EMsoftOO_SOURCE_DIR}/NamelistTemplates")
set(LIB_SEARCH_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if((EMsoftOO_ENABLE_HDF5_SUPPORT) AND (EMsoftOO_ENABLE_OpenCL_SUPPORT))

  GetHDF5LinkLibraries(EMSOFTOO)
  # The libraries are specifically laid out this way in order to ensure the link line
  # has the correct ordering. This seems to be more important on Windows than on macOS/Linux
  set(EXE_LINK_LIBRARIES  
          ${OpenCL_LIBRARY}
          clfortran
          EMsoftOOLib 
          EMOpenCLLib
          ${EMSOFTOO_hdf5LinkLibs}
  )

  Add_EMSoftOO_Executable(TARGET testMCOpenCLWrapper 
                        SOURCES ${APP_DIR}/testMCOpenCLWrapper.f90 
                        LINK_LIBRARIES ${EXE_LINK_LIBRARIES} 
                        TEMPLATE ${TMPLT_DIR}/EMMCOpenLC.template
                        SOLUTION_FOLDER EMsoftOOPublic/TestPrograms 
                        INSTALL_PROGRAM TRUE
                        INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR} ${EMOpenCLLib_BINARY_DIR}
                        )
endif()

if(EMsoftOO_ENABLE_HDF5_SUPPORT)

      GetHDF5LinkLibraries(EMSOFTOO)
      set(EXE_LINK_LIBRARIES ${EXE_LINK_LIBRARIES} ${EMSOFTOO_hdf5LinkLibs})

      Add_EMSoftOO_Executable(TARGET testEBSDmasterWrapper
                                   SOURCES ${APP_DIR}/testEBSDmasterWrapper.f90
                                   LINK_LIBRARIES EMsoftOOLib ${EMSOFTOO_hdf5LinkLibs}
                                   TEMPLATE ${TMPLT_DIR}/EMEBSDmaster.template
                                   SOLUTION_FOLDER EMsoftOOPublic/TestPrograms
                                   INSTALL_PROGRAM TRUE
                                   INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                                   )

      Add_EMSoftOO_Executable(TARGET tester
                                   SOURCES ${APP_DIR}/tester.f90
                                   LINK_LIBRARIES EMsoftOOLib ${EMSOFTOO_hdf5LinkLibs}
                                   TEMPLATE ${TMPLT_DIR}/EMEBSDmaster.template
                                   SOLUTION_FOLDER EMsoftOOPublic/TestPrograms
                                   INSTALL_PROGRAM TRUE
                                   INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                                   )

      Add_EMSoftOO_Executable(TARGET colortest
                                   SOURCES ${APP_DIR}/colortest.f90
                                   LINK_LIBRARIES EMsoftOOLib ${EMSOFTOO_hdf5LinkLibs}
                                   SOLUTION_FOLDER EMsoftOOPublic/TestPrograms
                                   INSTALL_PROGRAM TRUE
                                   INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                                   )

      Add_EMSoftOO_Executable(TARGET EMHSGtest
                                   SOURCES ${APP_DIR}/EMHSGtest.f90
                                   LINK_LIBRARIES EMsoftOOLib ${EMSOFTOO_hdf5LinkLibs}
                                   SOLUTION_FOLDER EMsoftOOPublic/TestPrograms
                                   INSTALL_PROGRAM TRUE
                                   INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                                   )



endif()
