

set_property(GLOBAL PROPERTY EMsoftOO_PACKAGE_DEST_PREFIX ".")
# -----------------------------------------------------------------------
#
# -----------------------------------------------------------------------

get_property(EMsoftOO_PACKAGE_DEST_PREFIX GLOBAL PROPERTY EMsoftOO_PACKAGE_DEST_PREFIX)

include("${EMsoftOO_SOURCE_DIR}/Source/EMsoftOO_Functions.cmake")

add_subdirectory(${PROJECT_SOURCE_DIR}/Source/EMsoftOOLib ${PROJECT_BINARY_DIR}/EMsoftOOLib)

option(EMsoftOO_ENABLE_PYTHON_SUPPORT "Build the C-interop shared library for Python bindings" OFF)
if( ${EMsoftOO_ENABLE_PYTHON_SUPPORT} )
  add_subdirectory(${PROJECT_SOURCE_DIR}/Source/EMsoftOOLib/c_interface ${PROJECT_BINARY_DIR}/EMsoftOO_c)
endif()

option(EMsoftOO_ENABLE_HDF5_SUPPORT "Enable HDF5 based I/O" ON)

option(EMsoftOO_ENABLE_OpenCL_SUPPORT "Enable OpenCL support" ON)

# Apple Metal GPU backend (see MetalMigrationPlan.md).  When ON, EMOpenCLLib is
# built with the Metal backend (mod_CLsupport_metal.f90 + metal-cpp shim) instead
# of the OpenCL one; the two expose the identical mod_CLsupport / OpenCL_T
# interface, so the program modules are unchanged.  Defaults ON on Apple (where
# Apple OpenCL is deprecated), OFF and forced off elsewhere.
if(APPLE)
  option(EMsoftOO_ENABLE_Metal_SUPPORT "Enable Apple Metal GPU backend (macOS)" ON)
  if(EMsoftOO_ENABLE_Metal_SUPPORT)
    # The *.metal -> *.metallib build step needs the Metal compiler
    # (`xcrun -sdk macosx metal`), which lives in full Xcode (and, on Xcode 16+,
    # its separately-downloaded Metal Toolchain component) -- it is NOT in a
    # Command-Line-Tools-only install.  If it is missing, fall back to the OpenCL
    # backend rather than failing the build, so a clean default build still works.
    execute_process(COMMAND xcrun -sdk macosx -f metal
                    RESULT_VARIABLE EMsoftOO_METAL_CC_RC OUTPUT_QUIET ERROR_QUIET)
    if(NOT EMsoftOO_METAL_CC_RC EQUAL 0)
      message(WARNING
        "The Metal compiler (`xcrun -sdk macosx metal`) was not found -- full "
        "Xcode (plus its Metal Toolchain component on Xcode 16+) is required to "
        "build the Metal backend. Falling back to the OpenCL GPU backend. Install "
        "the Metal toolchain and reconfigure with -DEMsoftOO_ENABLE_Metal_SUPPORT=ON "
        "to use Metal.")
      set(EMsoftOO_ENABLE_Metal_SUPPORT OFF CACHE BOOL "Enable Apple Metal GPU backend (macOS)" FORCE)
      # make sure a GPU backend still gets built (OpenCL) after the fallback
      if(NOT EMsoftOO_ENABLE_OpenCL_SUPPORT)
        set(EMsoftOO_ENABLE_OpenCL_SUPPORT ON CACHE BOOL "Enable OpenCL support" FORCE)
        message(STATUS "Re-enabling EMsoftOO_ENABLE_OpenCL_SUPPORT for the Metal->OpenCL fallback.")
      endif()
    endif()
  endif()
else()
  option(EMsoftOO_ENABLE_Metal_SUPPORT "Enable Apple Metal GPU backend (macOS)" OFF)
  set(EMsoftOO_ENABLE_Metal_SUPPORT OFF CACHE BOOL "Enable Apple Metal GPU backend (macOS)" FORCE)
endif()

if( ${EMsoftOO_ENABLE_OpenCL_SUPPORT} OR ${EMsoftOO_ENABLE_Metal_SUPPORT} )
  add_subdirectory(${PROJECT_SOURCE_DIR}/Source/EMOpenCLLib ${PROJECT_BINARY_DIR}/EMOpenCLLib)
endif()

set(MODALITY_DIRS
    CTEMbook
    Demag
    DictionaryIndexing
    GBs
    # EEC
    # OLIO
    OM
    #pyEMsoftOO
    SEM
    Shapes
    TEM
    SphInx
    QC
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
