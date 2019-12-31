#--////////////////////////////////////////////////////////////////////////////
# Copyright (c) 2009-2015 BlueQuartz Software, LLC
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright notice, this
# list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
#
# Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
# contributors may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# The code contained herein was partially funded by the followig contracts:
#    United States Air Force Prime Contract FA8650-07-D-5800
#    United States Air Force Prime Contract FA8650-10-D-5210
#    United States Prime Contract Navy N00173-07-C-2068
#--////////////////////////////////////////////////////////////////////////////
if(NOT DEFINED CMP_SOURCE_DIR)
  get_filename_component(CMP_SOURCE_DIR ${CMAKE_CURRENT_LIST_FILE} PATH CACHE)
endif()
set(CMP_CONFIGURED_FILES_SOURCE_DIR ${CMP_SOURCE_DIR}/ConfiguredFiles CACHE INTERNAL "")
set(CMP_CORE_TESTS_SOURCE_DIR ${CMP_SOURCE_DIR}/CoreTests CACHE INTERNAL "")
set(CMP_INSTALLATION_SUPPORT_SOURCE_DIR ${CMP_SOURCE_DIR}/InstallationSupport CACHE INTERNAL "")
set(CMP_MODULES_SOURCE_DIR ${CMP_SOURCE_DIR}/Modules CACHE INTERNAL "")
set(CMP_OSX_TOOLS_SOURCE_DIR ${CMP_SOURCE_DIR}/OSX_Tools CACHE INTERNAL "")
set(CMP_LINUX_TOOLS_SOURCE_DIR ${CMP_SOURCE_DIR}/Linux_Tools CACHE INTERNAL "")

# --------------------------------------------------------------------
# Over ride CMake's built in module directory by prepending cmp's module
# directory first
set(CMAKE_MODULE_PATH ${CMP_MODULES_SOURCE_DIR} ${CMAKE_MODULE_PATH})

INCLUDE(${CMP_SOURCE_DIR}/cmpCMakeMacros.cmake )

include( ${CMP_CORE_TESTS_SOURCE_DIR}/cmpConfigureChecks.cmake )
if(NOT DEFINED CMP_PROJECT_NAMESPACE)
    set(CMP_PROJECT_NAMESPACE "CMP")
endif()

if(NOT DEFINED CMP_HEADER_DIR)
    set(CMP_HEADER_DIR ${PROJECT_BINARY_DIR}/cmp)
endif()

if(NOT DEFINED CMP_CONFIGURATION_FILE_NAME)
    set(CMP_CONFIGURATION_FILE_NAME "cmpConfiguration.h")
endif()

if(NOT DEFINED CMP_TYPES_FILE_NAME)
    set(CMP_TYPES_FILE_NAME "cmpTypes.h")
endif()

if(NOT DEFINED CMP_VERSION_HEADER_FILE_NAME)
    set(CMP_VERSION_HEADER_FILE_NAME "cmpVersion.h")
endif()

if(NOT DEFINED CMP_VERSION_SOURCE_FILE_NAME)
    set(CMP_VERSION_HEADER_FILE_NAME "cmpVersion.cpp")
endif()

if(NOT DEFINED CMP_VERSION_MACRO_FILE_NAME)
    set(CMP_VERSION_MACRO_FILE_NAME "cmpVersionMacro.h")
endif()

if(NOT DEFINED CMP_TOP_HEADER_FILE)
    set(CMP_TOP_HEADER_FILE "")
endif()

if(NOT DEFINED CMP_PROJECT_NAME)
    set(CMP_PROJECT_NAME "CMP")
endif()

if(NOT DEFINED CMP_PLUGIN_LIST_FILE)
    set(CMP_PLUGIN_LIST_FILE ${PROJECT_BINARY_DIR}/plugins.txt)
endif()

if(NOT DEFINED CMP_PLUGIN_SEARCHDIR_FILE)
    set(CMP_PLUGIN_SEARCHDIR_FILE ${PROJECT_BINARY_DIR}/libsearchdirs.txt)
endif()

get_filename_component(CMP_CONFIGURATION_HEADER_GUARD ${CMP_CONFIGURATION_FILE_NAME} NAME_WE)
get_filename_component(CMP_TYPES_HEADER_GUARD ${CMP_TYPES_FILE_NAME} NAME_WE)
get_filename_component(CMP_VERSION_HEADER_GUARD ${CMP_VERSION_HEADER_FILE_NAME} NAME_WE)

# --------------------------------------------------------------------
# Generate our header files
# --------------------------------------------------------------------
cmpConfigureFileWithMD5Check(CONFIGURED_TEMPLATE_PATH ${CMP_CONFIGURED_FILES_SOURCE_DIR}/cmpConfiguration.h.in
                            GENERATED_FILE_PATH ${CMP_HEADER_DIR}/${CMP_CONFIGURATION_FILE_NAME} )
cmpConfigureFileWithMD5Check(CONFIGURED_TEMPLATE_PATH ${CMP_CONFIGURED_FILES_SOURCE_DIR}/cmpPrimitiveTypes.h.in
                             GENERATED_FILE_PATH ${CMP_HEADER_DIR}/${CMP_TYPES_FILE_NAME} )


# --------------------------------------------------------------------
# Generate a Header file with Compile Version variables
# --------------------------------------------------------------------
if( ${CMP_GENERATE_VERSION_STRING} )
        # Find Git executable
    Find_package(Git)

    set(PROJECT_EXPORT_MACRO ${CMP_PROJECT_NAMESPACE}_EXPORT)
    if( "${CMP_PROJECT_NAMESPACE}" STREQUAL "")
      set(PROJECT_EXPORT_MACRO "")
    endif()

    cmpRevisionString( GENERATED_HEADER_FILE_PATH "${CMP_VERSION_HEADER_FILE_NAME}"
                            GENERATED_SOURCE_FILE_PATH "${CMP_VERSION_SOURCE_FILE_NAME}"
                            GENERATED_MACRO_HEADER_PATH "${CMP_VERSION_MACRO_FILE_NAME}"
                            NAMESPACE "${CMP_PROJECT_NAMESPACE}"
                            PROJECT_NAME "${PROJECT_NAME}"
                            EXPORT_MACRO "${PROJECT_EXPORT_MACRO}")
   
endif()

cmp_IDE_GENERATED_PROPERTIES( "Generated"
              "${CMP_HEADER_DIR}/${CMP_CONFIGURATION_FILE_NAME}"
              "${CMP_HEADER_DIR}/${CMP_TYPES_FILE_NAME}"
              "${CMP_HEADER_DIR}/${CMP_VERSION_HEADER_FILE_NAME};${CMP_HEADER_DIR}/${CMP_VERSION_SOURCE_FILE_NAME}")

# --------------------------------------------------------------------
# Enable the use of plugins that will get generated as part of the project
# We are going to write the paths to the plugins into a file and then that
# file will be used as input to set an actual cmake variable and then
# passed to the bundle utilities cmake macro.
if(CMP_ENABLE_PLUGINS)

  file(WRITE ${CMP_PLUGIN_LIST_FILE} "")
  file(WRITE ${CMP_PLUGIN_SEARCHDIR_FILE} "${PROJECT_BINARY_DIR}/Bin/Plugins;")
  file(APPEND ${CMP_PLUGIN_SEARCHDIR_FILE} "${PROJECT_BINARY_DIR}/Bin;")

endif()



