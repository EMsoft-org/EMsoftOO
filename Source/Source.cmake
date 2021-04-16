

set_property(GLOBAL PROPERTY EMsoftOO_PACKAGE_DEST_PREFIX ".")
# -----------------------------------------------------------------------
#
# -----------------------------------------------------------------------

get_property(EMsoftOO_PACKAGE_DEST_PREFIX GLOBAL PROPERTY EMsoftOO_PACKAGE_DEST_PREFIX)

include("${EMsoftOO_SOURCE_DIR}/Source/EMsoftOO_Functions.cmake")

add_subdirectory(${PROJECT_SOURCE_DIR}/Source/EMsoftOOLib ${PROJECT_BINARY_DIR}/EMsoftOOLib)

option(EMsoftOO_ENABLE_HDF5_SUPPORT "Enable HDF5 based I/O" ON)

option(EMsoftOO_ENABLE_OpenCL_SUPPORT "Enable OpenCL support" ON)
if( ${EMsoftOO_ENABLE_OpenCL_SUPPORT} )
  add_subdirectory(${PROJECT_SOURCE_DIR}/Source/EMOpenCLLib ${PROJECT_BINARY_DIR}/EMOpenCLLib)
endif()

set(MODALITY_DIRS
    DictionaryIndexing
    # GBs
    # EEC
    # OLIO
    OM
    #pyEMsoftOO
    SEM
    # TEM
    # QC
    Utilities
    TestPrograms
    XRay
)
# -----------------------------------------------------------------------
# Establish which modalities are going to be compiled
# -----------------------------------------------------------------------
foreach(MODALITY ${MODALITY_DIRS})
  option(EMsoftOO_ENABLE_${MODALITY} "Build sources and programs related to ${MODALITY}" ON)
endforeach()


# -----------------------------------------------------------------------
# Add a wrapper lib thats uses the enabled modality options to compile itself
# -----------------------------------------------------------------------
# add_subdirectory(${PROJECT_SOURCE_DIR}/Source/EMsoftOOWrapperLib ${PROJECT_BINARY_DIR}/EMsoftOOWrapperLib)

# -----------------------------------------------------------------------
# Add the executables
# -----------------------------------------------------------------------
foreach(MODALITY ${MODALITY_DIRS})
  if( "${EMsoftOO_ENABLE_${MODALITY}}" STREQUAL "ON" )
    message(STATUS "EMsoftOO: Enabling public ${MODALITY} Modality")
    add_subdirectory( ${PROJECT_SOURCE_DIR}/Source/${MODALITY} ${PROJECT_BINARY_DIR}/${MODALITY})
  endif()
endforeach()



# -----------------------------------------------------------------------
# Does the developer want to compile the GUI for EMsoftOO?
# -----------------------------------------------------------------------
# if( EMsoftOO_ENABLE_EMsoftOOWorkbench )

#   INCLUDE (${EMsoftOO_SOURCE_DIR}/Support/cmp/cmpCMakeMacros.cmake )
#   # --------------------------------------------------------------------
#   # Find and Use the Qt5 Libraries
#   include(${EMsoftOO_SOURCE_DIR}/Support/cmp/ExtLib/Qt5Support.cmake)
#   set(EMsoftOOWorkbench_Qt5_Components Core Widgets Network Gui Concurrent Svg Xml OpenGL PrintSupport )
#   CMP_AddQt5Support( "${EMsoftOOWorkbench_Qt5_Components}"
#                     "FALSE"
#                     "${EMsoftOO_BINARY_DIR}"
#                     "EMsoftOOWorkbench")

#   include(${PROJECT_SOURCE_DIR}/Source/EMsoftOOWorkbench/SourceList.cmake)
# endif()
