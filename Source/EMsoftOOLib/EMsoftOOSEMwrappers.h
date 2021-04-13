#ifndef _emsoftoo_SEMwrappers_H_
#define _emsoftoo_SEMwrappers_H_


#ifdef __cplusplus
#include <cstddef>
#include <cstdint>
extern "C" {
#endif

/**
* @brief This is the typedef for a call back function that is
* used in the EMsoft library.
* @param size_t Unique integer that designates which C++ object
* did the call into EMsoft
* @param int
*/
typedef void (*ProgCallBackType)(size_t, int);

/**
* @brief This is the typedef for a call back function that is
* used in the EMsoft library.
* @param size_t Unique integer that designates which C++ object
* did the call into EMsoft
* @param int
* @param int
* @param float
*/
typedef void (*ProgCallBackType2)(size_t, int, int, float);

/**
* @brief This is the typedef for a call back function that is
* used in the EMsoft library.
* @param size_t Unique integer that designates which C++ object
* did the call into EMsoft
* @param int
* @param int
* @param int
* @param int
*/
typedef void (*ProgCallBackType3)(size_t, int, int, int, int);