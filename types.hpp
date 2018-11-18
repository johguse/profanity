#ifndef HPP_TYPES
#define HPP_TYPES

/* The structs declared in this file should have size/alignment hints
 * to ensure that their representation is identical to that in OpenCL.
 */
#include <CL/cl.h>

#define MP_NWORDS 8

typedef cl_uint mp_word;

typedef struct {
	mp_word d[MP_NWORDS];
} mp_number;

typedef struct {
    mp_number x;
    mp_number y;
} point;

typedef struct {
	cl_uint found;
	cl_uint foundId;
	cl_uchar foundHash[20];
} result;

#endif /* HPP_TYPES */