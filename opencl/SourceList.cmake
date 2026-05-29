
#---------------------------------------------------------------------
# Set some variables to shorten up the call to the function below
set(APP_DIR ${EMsoftOO_SOURCE_DIR}/opencl)

#---------------------------------------------------------------------
# Aggregate all the OpenCL files that are needed
set(EMSoft_CL_SRCS
  ${APP_DIR}/DictIndx.cl
  ${APP_DIR}/EMMC.cl
  ${APP_DIR}/EMMCfoil.cl
  ${APP_DIR}/EMMCxyz.cl
  ${APP_DIR}/MBmoduleOpenCL.cl
)


if(NOT EXISTS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl")
  file(MAKE_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl")
endif()

foreach(file ${EMSoft_CL_SRCS})
  file(COPY "${file}" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl/")
endforeach()


#---------------------------------------------------------------------
# This sets up the two variables install_dir and lib_install_dir
EMsoft_SetupInstallDirs()

#---------------------------------------------------------------------
# Create the Installation Rules
INSTALL(FILES ${EMSoft_CL_SRCS}
  COMPONENT Applications
  DESTINATION ${extra_install_dir}/opencl
)

#---------------------------------------------------------------------
# When the Metal backend is enabled, EMOpenCLLib builds the *.metallib files
# into ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl/ at build time (see
# Source/EMOpenCLLib/CMakeLists.txt).  Install them next to the .cl files so
# packaged/installed runs resolve them via OpenCLpathname.  FILES_MATCHING picks
# up whatever metallibs were actually built (no hard-coded list to keep in sync).
if(EMsoftOO_ENABLE_Metal_SUPPORT)
  INSTALL(DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/opencl/
    COMPONENT Applications
    DESTINATION ${extra_install_dir}/opencl
    FILES_MATCHING PATTERN "*.metallib"
  )
endif()
