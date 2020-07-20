
include ("@CMP_SOURCE_DIR@/cmpCMakeMacros.cmake")

#------------------------------------------------------------------------------
# Generate Fortran module file that contains all the String constants
set(FORTRAN_STRING_CONSTANTS "")

file(STRINGS "@EMsoftOOLib_SOURCE_DIR@/stringconstants.in" EMSOFTOO_STRING_CONSTANTS)
foreach(SC ${EMSOFTOO_STRING_CONSTANTS})
  string(FIND "${SC}" "!" COMMENT_POS)
  list(LENGTH SC LIST_LENGTH)
  string(FIND "${SC}" "[CATEGORY]" CATEGORY_POS)
  if(NOT CATEGORY_POS EQUAL -1)
    list(GET SC 1  VAR_VALUE)
    set(FORTRAN_STRING_CONSTANTS
      ${FORTRAN_STRING_CONSTANTS}
      "\n!------------------------------------------\n"
      "! ${VAR_VALUE}\n"
      "!------------------------------------------\n"
      )
  elseif(COMMENT_POS EQUAL -1 AND LIST_LENGTH EQUAL 2)
    list(GET SC 0  VAR_NAME)
    list(GET SC 1  VAR_VALUE)
    string(LENGTH ${VAR_VALUE} VAR_STR_LEN)
    math(EXPR VAR_STR_LEN ${VAR_STR_LEN}-2)
    set(FORTRAN_STRING_CONSTANTS
      ${FORTRAN_STRING_CONSTANTS}
      "character(${VAR_STR_LEN}), parameter     :: SC_${VAR_NAME} = ${VAR_VALUE}\n!DEC$ ATTRIBUTES DLLEXPORT :: SC_${VAR_NAME}\n"
      )
  endif()

endforeach()


string(REPLACE  ";" "" FORTRAN_STRING_CONSTANTS ${FORTRAN_STRING_CONSTANTS})
cmpConfigureFileWithMD5Check(CONFIGURED_TEMPLATE_PATH "@EMsoftOOLib_SOURCE_DIR@/stringconstants.in.f90"
                             GENERATED_FILE_PATH "@EMsoftOOLib_BINARY_DIR@/stringconstants.f90"
                             VERBOSE TRUE)
