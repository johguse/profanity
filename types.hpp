#ifndef HPP_TYPES
#define HPP_TYPES

/* The structs declared in this file should have size/alignment hints
 * to ensure that their representation is identical to that in OpenCL.
 */
#include <CL/cl.h>

#define BN_NWORDS 8

typedef cl_uint bn_word;

typedef struct {
	bn_word d[BN_NWORDS];
} bignum;

typedef struct {
    bignum x;
    bignum y;
} point;

typedef struct {
	cl_uint found;
	cl_uint foundId;
	cl_uchar foundScore;
	cl_uchar foundHash[20];
} result;

#endif /* HPP_TYPES */