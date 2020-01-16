include("${EMsoftOO_SOURCE_DIR}/Source/EMsoftOO_Functions.cmake")
set(APP_DIR "${EMsoftOO_SOURCE_DIR}/Source/Test")
set(TMPLT_DIR "${EMsoftOO_SOURCE_DIR}/NamelistTemplates")

set(EMsoftOOTestDir "${EMsoftOO_SOURCE_DIR}/Source/Test")


#------------------------------------------------------------------------------
# in order to add a new unit test the following needs to be performed:
# For a Unit Test called "foo":
# Create a fooTest.f90 file that implements the following function:
# module foo
#
# contains
#  subroutine fooExecuteTest(result) &
#  bind(c, name='fooExecuteTest') ! this routine is callable from a C/C++ program
#  !DEC$ ATTRIBUTES DLLEXPORT :: fooExecuteTest
#
# end subroutine fooExecuteTest
#
# end module foo
#
# Then add the proper cmake code below to instantiate the test:
# AddEMsoftOOUnitTest(TARGET fooTest
#                   SOURCES ${EMsoftOOTestDir}/fooTest.f90
#                   TEST_NAME fooTest
#                   LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
#                   SOLUTION_FOLDER EMsoftOOPublic/Test)
#
# The argument LINK_LIBRARIES needs to have all the proper libraries that the
# test will link against. EXE_LINK_LIBRARIES is just a CMake variable that
# the EMsoftOO project uses internally to list the most common 3rd part and
# internal libraries that are used. Your test may have additional libraries
# such as JSON Fortran or BLAS/LAPACK for example that are also needed.

#  /NODEFAULTLIB:msvcrt.lib /NODEFAULTLIB:msvcmrt.lib /NODEFAULTLIB:msvcurt.lib /NODEFAULTLIB:msvcrtd.lib

if(EMsoftOO_ENABLE_HDF5_SUPPORT)
  
  GetHDF5LinkLibraries(EMSOFTOO)
  set(EXE_LINK_LIBRARIES ${EMSOFTOO_hdf5LinkLibs} EMsoftOOLib )

  set(EMsoftOO_UNIT_TESTS 
    MODQuaternions
    HDFtextfile
    HDFinteger
    HDFchar
    HDFfloat
    HDFdouble
    HDFhyperslab
    MODsymmetry
    Image
    # MODcrystal
    # MODMuellerCalculus
    # MODRotations
  )

  foreach(test ${EMsoftOO_UNIT_TESTS})
    message(STATUS "Adding Unit test: ${test}" )
    AddEMsoftOOUnitTest(TARGET ${test}Test
                      SOURCES ${EMsoftOOTestDir}/${test}Test.f90
                      TEST_NAME ${test}
                      LINK_LIBRARIES ${EXE_LINK_LIBRARIES}
                      SOLUTION_FOLDER EMsoftOOPublic/Test
                      INCLUDE_DIRS ${EMsoftOO_SOURCE_DIR}/Source/Test ${EMsoftOOLib_BINARY_DIR})
    if(MSVC)
        target_compile_options( ${test}Test PRIVATE "/MT$<$<CONFIG:Debug>:d>")
        target_link_options(${test}Test PRIVATE "/NODEFAULTLIB:MSVCRT")
    endif()
  endforeach()

endif()

#------------------------------------------------------------------------------
# Create CTest entries for a basic indexing workflow
#------------------------------------------------------------------------------
# Use multiple CPU cores to Test
include(ProcessorCount)
if(NOT DEFINED EMsoftOO_TEST_CPU_COUNT) 
  ProcessorCount(EMsoftOO_TEST_CPU_COUNT)
endif()
# message(STATUS "Using ${EMsoftOO_TEST_CPU_COUNT} CPU Cores for tests")

set(EXE_EXT "")
if(WIN32)
  set(EXE_EXT ".exe")
endif()

#------------------------------------------------------------------------------
# This bit ensures that the EMsoftOOData_DIR is set correctly or errors out.
get_filename_component(EMSOFTOO_PARENT ${EMsoftOO_SOURCE_DIR} DIRECTORY CACHE)
if("${EMsoftOOData_DIR}" STREQUAL "")
  if(EXISTS "${EMSOFTOO_PARENT}/EMsoftData")
    set(EMsoftOOData_DIR ${EMSOFTOO_PARENT}/EMsoftData)
  else()
    message(WARNING "EMsoftData_DIR not found or empty folder, cloning automatically from https://github.com/emsoft-org/EMsoftData")
    message(FATAL_ERROR "${EMSOFTOO_PARENT}/EMsoftData")
    include(ExternalProject)
    ExternalProject_Add("EMsoftOODataClone"
      GIT_REPOSITORY https://github.com/emsoft-org/EMsoftData
      GIT_TAG "develop" # get the develop branch
      DOWNLOAD_DIR "${EMSOFTOO_PARENT}/EMsoftData"
      GIT_SHALLOW TRUE # we don't need version history
      GIT_PROGRESS TRUE # there may be some large datafiles
      # CONFIGURE_COMMAND "" # nothing needs to be configured
      # BUILD_COMMAND "" # nothing needs to be built
      # INSTALL_COMMAND "" # nothing needs to be installed
      # TEST_COMMAND "" # nothing needs to be tested
    )

    # message(FATAL_ERROR "EMsoftOOData CMake Variable must point to the path to the EMsoftOOData repository.\
    # We looked in the same directory as the 'EMsoftOO' directory for a directory called EMsoftOOData but did not find it.\
    # You can clone it from https://github.com/emsoft-org/EMsoftOOData using the 'develop' branch\
    # or you can set the EMsoftOOData_DIR variable to point to the location of the EMsoftOOData directory.\    
    # ")  
  endif()
endif()

#------------------------------------------------------------------------------
#  Function: create_workflow_test
#  brief: This creates a ctest for a specific EBSD DI Workflow
#
# function(create_workflow_test)
#   set(options )
#   set(oneValueArgs EXE_NAME INDEX OUTPUT_DIR EXTRA_NML)
#   set(multiValueArgs )
#   cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

#   if(NOT EXISTS ${EMsoftOOData_DIR}/DItutorial/Ni/${TEST_EXE_NAME}.nml)
#     message(FATAL_ERROR "NML File does not exist: ${EMsoftOOData_DIR}/DItutorial/Ni/${TEST_EXE_NAME}.nml\
#     This file is part of the EMsoftOOData repository. Please ensure that the EMsoftOOData repository is\
#     up to date via a 'git pull origin develop' command. The current value of the EMsoftOOData_DIR is \
#     ${EMsoftOOData_DIR}")
#   endif()

#   configure_file(${EMsoftOOData_DIR}/DItutorial/Ni/${TEST_EXE_NAME}.nml
#                   ${TEST_OUTPUT_DIR}/${TEST_EXE_NAME}.nml
#   )

#   if(NOT "${TEST_EXTRA_NML}" STREQUAL "")
#     if(NOT EXISTS ${EMsoftOOData_DIR}/DItutorial/Ni/${TEST_EXTRA_NML})
#       message(FATAL_ERROR "NML File does not exist: ${EMsoftOOData_DIR}/DItutorial/Ni/${TEST_EXTRA_NML}\
#       This file is part of the EMsoftOOData repository. Please ensure that the EMsoftOOData repository is\
#       up to date via a 'git pull origin develop' command. The current value of the EMsoftOOData_DIR is \
#       ${EMsoftOOData_DIR}")
#     endif()
  
#     configure_file(${EMsoftOOData_DIR}/DItutorial/Ni/${TEST_EXTRA_NML}
#                   ${EMsoftOO_NML_TEST_DIR}/${TEST_EXTRA_NML}
#     )
#   endif()
#   message(STATUS "Adding Unit Test: EMsoftOO_${TEST_INDEX}_${TEST_EXE_NAME}")
#   add_test(NAME EMsoftOO_${TEST_INDEX}_${TEST_EXE_NAME}
#           COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TEST_EXE_NAME}${EXE_EXT}" "${TEST_OUTPUT_DIR}/${TEST_EXE_NAME}.nml" 
#           WORKING_DIRECTORY ${EMsoftOO_NML_TEST_DIR})

# endfunction()

#------------------------------------------------------------------------------
# If the developer did not set these values then we need something.
if("${EMsoftOO_TEST_GPU_PLATFORM_ID}" STREQUAL "")
  set(EMsoftOO_TEST_GPU_PLATFORM_ID 1)
endif()

if("${EMsoftOO_TEST_GPU_DEVICE_ID}" STREQUAL "")
  set(EMsoftOO_TEST_GPU_DEVICE_ID 1)
endif()



configure_file("${EMsoftOO_SOURCE_DIR}/Source/Test/TestFileLocations.in.h"
  ${EMsoftOO_BINARY_DIR}/Test/TestFileLocations.h @ONLY
  )

# add_executable(DIWorkflowTestPrep DIWorkflowTestPrep.cpp ${EMsoftOO_BINARY_DIR}/Test/TestFileLocations.h)
# target_compile_features(DIWorkflowTestPrep PUBLIC cxx_std_11)
# target_include_directories(DIWorkflowTestPrep PUBLIC ${EMsoftOO_BINARY_DIR}/Test)

# #------------------------------------------------------------------------------
# # 1: The Ni.xtal file should have been copied into the proper spot via some other external means
# # We are just going to have to trust that it is there. We are numbering the tests
# # so that CTest runs them in the correct order.
# set(test_index  "1")
# set(EMsoftOO_NML_TEST_DIR "${EMsoftOO_BINARY_DIR}/DItutorial/Ni")
# file(MAKE_DIRECTORY ${EMsoftOO_NML_TEST_DIR})
# message(STATUS "Adding Unit Test: EMsoftOO_${test_index}_DIWorkflowTestPrep")
# add_test(NAME EMsoftOO_${test_index}_DIWorkflowTestPrep
#         COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/DIWorkflowTestPrep${EXE_EXT}
#         WORKING_DIRECTORY ${EMsoftOO_NML_TEST_DIR}
#         )

# #------------------------------------------------------------------------------
# # 2: Setup the EMMCOpenCL executable
# math(EXPR test_index "${test_index} + 1")
# create_workflow_test( EXE_NAME EMMCOpenCL
#                       INDEX ${test_index}
#                       OUTPUT_DIR ${EMsoftOO_NML_TEST_DIR}
#                     )


# #------------------------------------------------------------------------------
# # 3: Setup the EMEBSDmaster executable
# math(EXPR test_index "${test_index} + 1")
# create_workflow_test( EXE_NAME EMEBSDmaster
#                       INDEX ${test_index}
#                       OUTPUT_DIR ${EMsoftOO_NML_TEST_DIR}
#                       EXTRA_NML BetheParameters.nml
#                     )
# math(EXPR test_index "${test_index} + 1")
# create_workflow_test( EXE_NAME EMEBSD
#                       INDEX ${test_index}
#                       OUTPUT_DIR ${EMsoftOO_NML_TEST_DIR}
#                     )

# #------------------------------------------------------------------------------
# # 4: Setup the ADP executable
# math(EXPR test_index "${test_index} + 1")
# create_workflow_test( EXE_NAME EMgetADP
#                       INDEX ${test_index}
#                       OUTPUT_DIR ${EMsoftOO_NML_TEST_DIR}
#                     )

# #------------------------------------------------------------------------------
# # 5: Setup the EMEBSDDIpreview executable
# math(EXPR test_index "${test_index} + 1")
# create_workflow_test( EXE_NAME EMEBSDDIpreview
#                       INDEX ${test_index}
#                       OUTPUT_DIR ${EMsoftOO_NML_TEST_DIR}
#                     )


# #------------------------------------------------------------------------------
# # 6: Setup the EMEBSDDIpreview executable
# math(EXPR test_index "${test_index} + 1")
# create_workflow_test( EXE_NAME EMEBSDDI
#                       INDEX ${test_index}
#                       OUTPUT_DIR ${EMsoftOO_NML_TEST_DIR}
#                     )

# #------------------------------------------------------------------------------
# # 7: Setup the EMFitOrientation executable
# math(EXPR test_index "${test_index} + 1")
# create_workflow_test( EXE_NAME EMFitOrientation
#                       INDEX ${test_index}
#                       OUTPUT_DIR ${EMsoftOO_NML_TEST_DIR}
#                     )

