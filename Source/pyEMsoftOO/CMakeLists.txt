

# for now (Release 5.0), the python wrapper build is only done for Mac OS X or Linux
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/pyEMsoftOO")
  set(USERHOMEDIR "$ENV{HOME}")
  get_filename_component(EMsoftOO_BUILD_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} DIRECTORY)
  #-------------------------------------------------------------------------------
  # add run_pyEMsoftOO.sh to Build folder with properly substituted folder locations
  Message(STATUS "Configuring run_pyEMsoftOO.sh script in ${EMsoftOO_BUILD_DIR} folder")
  cmpConfigureFileWithMD5Check( GENERATED_FILE_PATH        ${EMsoftOO_BUILD_DIR}/run_pyEMsoftOO.sh
                                CONFIGURED_TEMPLATE_PATH   ${APP_DIR}/run_pyEMsoftOO.sh.in  )
  cmpConfigureFileWithMD5Check( GENERATED_FILE_PATH        ${EMsoftOO_BUILD_DIR}/run_docs.sh
                                CONFIGURED_TEMPLATE_PATH   ${APP_DIR}/run_docs.sh.in  )
endif()

if (Fortran_COMPILER_NAME MATCHES "ifort.*")
  if(WIN32)
    #Message(STATUS "Python wrapper script is currently not available on Windows")
      set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/pyEMsoftOO")
    #set(USERHOMEDIR "$ENV{HOME}")
    get_filename_component(EMsoftOO_BUILD_DIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} DIRECTORY)
    #-------------------------------------------------------------------------------
    # add run_pyEMsoftOO.sh to Build folder with properly substituted folder locations
    Message(STATUS "Configuring run_pyEMsoftOO.sh script in ${EMsoftOO_BUILD_DIR} folder")
    cmpConfigureFileWithMD5Check( GENERATED_FILE_PATH        ${EMsoftOO_BUILD_DIR}/run_pyEMsoftOO.sh
                                  CONFIGURED_TEMPLATE_PATH   ${APP_DIR}/run_pyEMsoftOO.sh.in  )
    #cmpConfigureFileWithMD5Check( GENERATED_FILE_PATH        ${EMsoftOO_BUILD_DIR}/run_docs.sh
    #                              CONFIGURED_TEMPLATE_PATH   ${APP_DIR}/run_docs.sh.in  )
  endif()
endif()

