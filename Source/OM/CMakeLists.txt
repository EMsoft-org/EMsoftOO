set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/OM")
set(TMPLT_DIR "${EMsoftOO_SOURCE_DIR}/NamelistTemplates")
set(LIB_SEARCH_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if(EMsoftOO_ENABLE_HDF5_SUPPORT)

  GetHDF5LinkLibraries(EMSOFTOO)
  set(EXE_LINK_LIBRARIES ${EXE_LINK_LIBRARIES} ${EMSOFTOO_hdf5LinkLibs} )

  Add_EMsoftOO_Executable(TARGET EMCPLMmaster
                              SOURCES ${APP_DIR}/EMCPLMmaster.f90 
                              LINK_LIBRARIES jsonfortran ${EXE_LINK_LIBRARIES}  
                              TEMPLATE ${TMPLT_DIR}/EMCPLMmaster.template 
                              SOLUTION_FOLDER EMsoftOOPublic/SEM/OM
                              INSTALL_PROGRAM TRUE
                              INCLUDE_DIRS ${EMsoftOOLib_BINARY_DIR}
                              )


endif()