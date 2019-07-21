/* ------------------------------------------------------------------------ */
/* Multiprecision functions                                                 */
/* ------------------------------------------------------------------------ */
#define MP_WORDS 8
#define MP_BITS 32
#define bswap32(n) (rotate(n & 0x00FF00FF, 24U)|(rotate(n, 8U) & 0x00FF00FF))

typedef uint mp_word;
typedef struct {
	mp_word d[MP_WORDS];
} mp_number;

__constant mp_number mod = { { 0xfffffc2f, 0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff} };
__constant mp_word mPrime = 0xd2253531;

mp_word mp_sub( mp_number * const r, const mp_number * const a, const mp_number * const b ) {
	mp_word t, c = 0;
	
	for( mp_word i = 0; i < MP_WORDS; ++i ) {
		t = a->d[i] - b->d[i] - c;
		c = t > a->d[i] ? 1 : (t == a->d[i] ? c : 0);
		
		r->d[i] = t;
	}
	
	return c;
}

mp_word mp_sub_mod( mp_number * const r ) {
	mp_word t, c = 0;
	
	for( mp_word i = 0; i < MP_WORDS; ++i ) {
		t = r->d[i] - mod.d[i] - c;
		c = t > r->d[i] ? 1 : (t == r->d[i] ? c : 0);
		
		r->d[i] = t;
	}
	
	return c;
}

void mp_mod_sub( mp_number * const r, const mp_number * const a, const mp_number * const b ) {
	mp_word i, t, c = 0;
	
	for( i = 0; i < MP_WORDS; ++i ) {
		t = a->d[i] - b->d[i] - c;
		c = t < a->d[i] ? 0 : (t == a->d[i] ? c : 1);
		
		r->d[i] = t;
	}
	
	if(c) {
		c = 0;
		for( i = 0; i < MP_WORDS; ++i ) {
			r->d[i] += mod.d[i] + c;
			c = r->d[i] < mod.d[i] ? 1 : (r->d[i] == mod.d[i] ? c : 0);
		}
	}
}

void mp_mod_sub_gx( mp_number * const r, const mp_number * const a) {
    // gx = {0x487e2097, 0xd7362e5a, 0x29bc66db, 0x231e2953, 0x33fd129c, 0x979f48c0, 0xe9089f48, 0x9981e643}
	//point g = { {  }, { {0xd3dbabe2, 0xb15ea6d2, 0x1f1dc64d, 0x8dfc5d5d, 0xac19c136, 0x70b6b59a, 0xd4a582d6, 0xcf3f851f} } };

	mp_word i, t, c = 0;
	
	t = a->d[0] - 0x487e2097; c = t < a->d[0] ? 0 : (t == a->d[0] ? c : 1); r->d[0] = t;
	t = a->d[1] - 0xd7362e5a - c; c = t < a->d[1] ? 0 : (t == a->d[1] ? c : 1); r->d[1] = t;
	t = a->d[2] - 0x29bc66db - c; c = t < a->d[2] ? 0 : (t == a->d[2] ? c : 1); r->d[2] = t;
	t = a->d[3] - 0x231e2953 - c; c = t < a->d[3] ? 0 : (t == a->d[3] ? c : 1); r->d[3] = t;
	t = a->d[4] - 0x33fd129c - c; c = t < a->d[4] ? 0 : (t == a->d[4] ? c : 1); r->d[4] = t;
	t = a->d[5] - 0x979f48c0 - c; c = t < a->d[5] ? 0 : (t == a->d[5] ? c : 1); r->d[5] = t;
	t = a->d[6] - 0xe9089f48 - c; c = t < a->d[6] ? 0 : (t == a->d[6] ? c : 1); r->d[6] = t;
	t = a->d[7] - 0x9981e643 - c; c = t < a->d[7] ? 0 : (t == a->d[7] ? c : 1); r->d[7] = t;
	
	if(c) {
		c = 0;
		for( i = 0; i < MP_WORDS; ++i ) {
			r->d[i] += mod.d[i] + c;
			c = r->d[i] < mod.d[i] ? 1 : (r->d[i] == mod.d[i] ? c : 0);
		}
	}
}

void mp_mod_sub_gy( mp_number * const r, const mp_number * const a) {
	mp_word i, t, c = 0;
	
	t = a->d[0] - 0xd3dbabe2; c = t < a->d[0] ? 0 : (t == a->d[0] ? c : 1); r->d[0] = t;
	t = a->d[1] - 0xb15ea6d2 - c; c = t < a->d[1] ? 0 : (t == a->d[1] ? c : 1); r->d[1] = t;
	t = a->d[2] - 0x1f1dc64d - c; c = t < a->d[2] ? 0 : (t == a->d[2] ? c : 1); r->d[2] = t;
	t = a->d[3] - 0x8dfc5d5d - c; c = t < a->d[3] ? 0 : (t == a->d[3] ? c : 1); r->d[3] = t;
	t = a->d[4] - 0xac19c136 - c; c = t < a->d[4] ? 0 : (t == a->d[4] ? c : 1); r->d[4] = t;
	t = a->d[5] - 0x70b6b59a - c; c = t < a->d[5] ? 0 : (t == a->d[5] ? c : 1); r->d[5] = t;
	t = a->d[6] - 0xd4a582d6 - c; c = t < a->d[6] ? 0 : (t == a->d[6] ? c : 1); r->d[6] = t;
	t = a->d[7] - 0xcf3f851f - c; c = t < a->d[7] ? 0 : (t == a->d[7] ? c : 1); r->d[7] = t;
	
	if(c) {
		c = 0;
		for( i = 0; i < MP_WORDS; ++i ) {
			r->d[i] += mod.d[i] + c;
			c = r->d[i] < mod.d[i] ? 1 : (r->d[i] == mod.d[i] ? c : 0);
		}
	}
}

mp_word mp_add( mp_number * const r, const mp_number * const a ) {
	mp_word c = 0;
	
	for( mp_word i = 0; i < MP_WORDS; ++i ) {
		r->d[i] += a->d[i] + c;
		c = r->d[i] < a->d[i] ? 1 : (r->d[i] == a->d[i] ? c : 0);
	}
	
	return c;
}

mp_word mp_add_mod( mp_number * const r ) {
	mp_word c = 0;
	
	for( mp_word i = 0; i < MP_WORDS; ++i ) {
		r->d[i] += mod.d[i] + c;
		c = r->d[i] < mod.d[i] ? 1 : (r->d[i] == mod.d[i] ? c : 0);
	}
	
	return c;
}

mp_word mp_add_more( mp_number * const r, mp_word * const extraR, const mp_number * const a, const mp_word * const extraA ) {
	const mp_word c = mp_add(r, a);
	*extraR += *extraA + c;
	return *extraR < *extraA ? 1 : (*extraR == *extraA ? c : 0);
}

mp_word mp_add_extra( mp_number * const r, mp_number * const a, mp_word * const extra ) {
	const mp_word c = mp_add(r, a);
	*extra += c;
	return *extra < c;
}

mp_word mp_gte( const mp_number * const a, const mp_number * const b ) {
	mp_word l = 0, g = 0;

	for( mp_word i = 0; i < MP_WORDS; ++i ) {
		if (a->d[i] < b->d[i]) l |= (1 << i);
		if (a->d[i] > b->d[i]) g |= (1 << i);
	}

	return g >= l;
}

void mp_shr_extra( mp_number * const r, mp_word * const e ) {
	r->d[0] = (r->d[1] << 31) | (r->d[0] >> 1);
	r->d[1] = (r->d[2] << 31) | (r->d[1] >> 1);
	r->d[2] = (r->d[3] << 31) | (r->d[2] >> 1);
	r->d[3] = (r->d[4] << 31) | (r->d[3] >> 1);
	r->d[4] = (r->d[5] << 31) | (r->d[4] >> 1);
	r->d[5] = (r->d[6] << 31) | (r->d[5] >> 1);
	r->d[6] = (r->d[7] << 31) | (r->d[6] >> 1);
	r->d[7] = (*e << 31)      | (r->d[7] >> 1);
	*e >>= 1;
}

void mp_shr( mp_number * const r ) {
	r->d[0] = (r->d[1] << 31) | (r->d[0] >> 1);
	r->d[1] = (r->d[2] << 31) | (r->d[1] >> 1);
	r->d[2] = (r->d[3] << 31) | (r->d[2] >> 1);
	r->d[3] = (r->d[4] << 31) | (r->d[3] >> 1);
	r->d[4] = (r->d[5] << 31) | (r->d[4] >> 1);
	r->d[5] = (r->d[6] << 31) | (r->d[5] >> 1);
	r->d[6] = (r->d[7] << 31) | (r->d[6] >> 1);
	r->d[7] >>= 1;
}

mp_word mp_mul_word( mp_number * const r, const mp_number * const a, const mp_word w, mp_word * const extra) {
	mp_word c = 0;
	for( mp_word i = 0; i < MP_WORDS; ++i ) {
		r->d[i] = a->d[i] * w + c;
		c = r->d[i] < c ? mul_hi(a->d[i], w) + 1 : mul_hi(a->d[i], w);
	}
	
	*extra += c;
	return *extra < c;
}

void mp_mul_mont( mp_number * const r, const mp_number * const a, const mp_number * const b) {
	mp_number mod_priv = { { 0xfffffc2f, 0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff} };
	mp_number A = { { 0, 0, 0, 0, 0, 0, 0, 0} };
	mp_number tmpNumber;
	
	mp_word extraWord = 0;
	mp_word overflow = 0;
	
	for( mp_word i = 0; i < MP_WORDS; ++i ) {
		/* Overflow used as temporary variable before being reset below. */
		overflow = (A.d[0] + a->d[i] * b->d[0]) * mPrime; // % b, where b = 2**MP_BITS
		overflow = mp_mul_word( &tmpNumber, &mod_priv, overflow, &extraWord);
		overflow += mp_add_extra( &A, &tmpNumber, &extraWord );
		overflow += mp_mul_word( &tmpNumber, b, a->d[i], &extraWord);
		overflow += mp_add_extra( &A, &tmpNumber, &extraWord );
	
		A.d[0] = A.d[1];
		A.d[1] = A.d[2];
		A.d[2] = A.d[3];
		A.d[3] = A.d[4];
		A.d[4] = A.d[5];
		A.d[5] = A.d[6];
		A.d[6] = A.d[7];
		A.d[7] = extraWord;
		extraWord = overflow;
	}
	
	if( extraWord ) { /* Ignore where N <= A < 2 ** 256 */
		mp_sub_mod(&A);
	}
	
	*r = A;
}

/* This is an optimization of the above. Moving out of Montgomery form can be done
 * by doing a Montgomery multiplication with 1. */
void mp_mul_mont_one( mp_number * const r, const mp_number * const a) {
	mp_number mod_priv = { { 0xfffffc2f, 0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff} };
	mp_number A = { { 0, 0, 0, 0, 0, 0, 0, 0} };
	mp_number tmpNumber;
	
	mp_word extraWord = 0;
	mp_word overflow = 0;
	mp_word c = 0;
	
	for( mp_word i = 0; i < MP_WORDS; ++i ) {
		/* Overflow used as temporary variable before being reset below. */
		overflow = (A.d[0] + a->d[i]) * mPrime; // % b, where b = 2**MP_BITS
		overflow = mp_mul_word( &tmpNumber, &mod_priv, overflow, &extraWord);
		overflow += mp_add_extra( &A, &tmpNumber, &extraWord );
		
		c = (A.d[0] + a->d[i]) < A.d[0];
		A.d[0] = A.d[1] + c; c = !A.d[0];
		A.d[1] = A.d[2] + c; c = !A.d[1];
		A.d[2] = A.d[3] + c; c = !A.d[2];
		A.d[3] = A.d[4] + c; c = !A.d[3];
		A.d[4] = A.d[5] + c; c = !A.d[4];
		A.d[5] = A.d[6] + c; c = !A.d[5];
		A.d[6] = A.d[7] + c; c = !A.d[6];
		A.d[7] = extraWord + c;
		extraWord = overflow;
	}
	
	if( extraWord ) { /* Ignore where N <= A < 2 ** 256 */
		mp_sub_mod(&A);
	}
	
	*r = A;
}

void mp_mod_inverse( mp_number * const r ) {
	mp_number A = { { 1 } };
	mp_number C = { { 0 } };
	mp_number v = mod;
	
	mp_word extraA = 0;
	mp_word extraC = 0;
	
	while( r->d[0] || r->d[1] || r->d[2] || r->d[3] || r->d[4] || r->d[5] || r->d[6] || r->d[7] ) {
		while( !(r->d[0] & 1) ) {
			mp_shr(r);
			if( A.d[0] & 1 ) {
				extraA += mp_add_mod(&A);
			}
				
			mp_shr_extra(&A, &extraA);
		}
		
		while( !(v.d[0] & 1) ) {
			mp_shr(&v);
			if( C.d[0] & 1 ) {
				extraC += mp_add_mod(&C);
			}
				
			mp_shr_extra(&C, &extraC);
		}
		
		if( mp_gte(r, &v) ) {
			mp_sub( r, r, &v );
			mp_add_more( &A, &extraA, &C, &extraC );
		} else {
			mp_sub(&v, &v, r);
			mp_add_more( &C, &extraC, &A, &extraA );
		}
	}
	
	while( extraC ) {
		extraC -= mp_sub_mod(&C);
	}
	




	v = mod;
	mp_sub(r, &v, &C);
}

/* ------------------------------------------------------------------------ */
/* Elliptic point and addition (with caveats).                              */
/* ------------------------------------------------------------------------ */
typedef struct {
    mp_number x;
    mp_number y;
} point;

// OpenCL crashes when trying to initialize local variables using the below
//__constant point generator = { { {0x487e2097, 0xd7362e5a, 0x29bc66db, 0x231e2953, 0x33fd129c, 0x979f48c0, 0xe9089f48, 0x9981e643} }, { {0xd3dbabe2, 0xb15ea6d2, 0x1f1dc64d, 0x8dfc5d5d, 0xac19c136, 0x70b6b59a, 0xd4a582d6, 0xcf3f851f} } };

// Does not handle points sharing X coordinate, this is a deliberate design choice.
void point_add(point * const p, point * const o) {
	mp_number mont_rrr = { { 0x3795f671, 0x002bb1e3, 0x00000b73, 0x1, 0, 0, 0, 0} };
	
	mp_number tmp;
	mp_number newX;
	mp_number newY;

	mp_mod_sub( &tmp, &o->x, &p->x );

	mp_mod_inverse( &tmp );
	mp_mul_mont( &tmp, &tmp, &mont_rrr);

	mp_mod_sub( &newX, &o->y, &p->y );
	mp_mul_mont( &tmp, &tmp, &newX );

	mp_mul_mont( &newX, &tmp, &tmp );
	mp_mod_sub( &newX, &newX, &p->x );
	mp_mod_sub( &newX, &newX, &o->x );

	mp_mod_sub( &newY, &p->x, &newX );
	mp_mul_mont( &newY, &newY, &tmp );
	mp_mod_sub( &newY, &newY, &p->y );

	p->x = newX;
	p->y = newY;
}

/* ------------------------------------------------------------------------ */
/* Profanity.                                                               */
/* ------------------------------------------------------------------------ */
typedef struct {
	uint found;
	uint foundId;
	uchar foundHash[20];
} result;

void profanity_begin_seed(__global const point * const precomp, point * const p, bool * const pIsFirst, const size_t precompOffset, const ulong seed) {
	point o;

	for( uchar i = 0; i < 8; ++i ) {
		const uchar shift = i * 8;
		const uchar byte = (seed >> shift) & 0xFF;

		if(byte) {
			o = precomp[precompOffset + i * 255 + byte - 1];
			if( *pIsFirst ) {
				*p = o;
				*pIsFirst = false;
			} else {
				point_add(p, &o);
			}
		}
	}
}

__kernel void profanity_begin(__global const point * const precomp, __global point * const pPoints, __global result * const pResult, const ulong4 seed) {
	const size_t id = get_global_id(0);
	point p;
	bool bIsFirst = true;

	profanity_begin_seed(precomp, &p, &bIsFirst, 8 * 255 * 0, seed.x);
	profanity_begin_seed(precomp, &p, &bIsFirst, 8 * 255 * 1, seed.y);
	profanity_begin_seed(precomp, &p, &bIsFirst, 8 * 255 * 2, seed.z);
	profanity_begin_seed(precomp, &p, &bIsFirst, 8 * 255 * 3, seed.w + id);

	pPoints[id] = p;

	for( uchar i = 0; i < PROFANITY_MAX_SCORE + 1; ++i ) {
		pResult[i].found = 0;
	}
}

__kernel void profanity_inverse_multiple(__global point * const pPoints, __global mp_number * const pInverse) {

	const size_t id = get_global_id(0) * PROFANITY_INVERSE_SIZE;
	
	mp_number copy1, copy2, copy3;
	mp_number buffer[PROFANITY_INVERSE_SIZE];
	mp_number mont_rrr = { { 0x3795f671, 0x002bb1e3, 0x00000b73, 0x1, 0, 0, 0, 0 } };

	buffer[0] = pPoints[id].x;
	mp_mod_sub_gx( &buffer[0], &buffer[0]);
	pInverse[0] = buffer[0];

	for( uint i = 1; i < PROFANITY_INVERSE_SIZE; ++i ) {
		buffer[i] = pPoints[id + i].x;
		mp_mod_sub_gx(&buffer[i], &buffer[i]);
		pInverse[id + i] = buffer[i];
		mp_mul_mont(&buffer[i], &buffer[i - 1], &buffer[i]);
	}

	// mp_mod_inverse(aR) -> (aR)^-1 
	// mp_mul_mont(x,y) -> x * y * R^-1
	// mp_mul_mont( (aR)^-1, R^3 ) -> (aR)^-1 * R^3 * R^-1 = a^-1 * R
	copy1 = buffer[PROFANITY_INVERSE_SIZE - 1];
	mp_mod_inverse( &copy1 );
	mp_mul_mont(&copy1, &copy1, &mont_rrr); 

	for( uint i = PROFANITY_INVERSE_SIZE - 1; i > 0; --i ) {
		copy2 = pInverse[id + i];

		mp_mul_mont(&copy3, &buffer[i - 1], &copy1);
		mp_mul_mont(&copy1, &copy2, &copy1);
		pInverse[id + i] = copy3;
	}

	pInverse[id] = copy1;
}

/*
// Unrolled version of the inversion algorithm where PROFANITY_INVERSE_SIZE = 64.
// This one gave horribly performance on my GTX 1070 and massively increased build time.
// On an RX480 it gave worse performance, 64MH/s vs 74MH/s.
// I'll leave it as a multi-line comment here for possible future experimentation on other platforms.

#define INVERSE_BEGIN(B) B = pPoints[id].x; mp_mod_sub_gx(&B, &B); pInverse[0] = B;
#define INVERSE(B, BP, offset) B = pPoints[id+offset].x; mp_mod_sub_gx(&B, &B); pInverse[id+offset] = B; mp_mul_mont(&B, &BP, &B);
#define INVERSE_REV(B, offset) copy2 = pInverse[id + offset]; mp_mul_mont(&copy3, &B, &copy1); mp_mul_mont(&copy1, &copy2, &copy1); pInverse[id + offset] = copy3;
#define INVERSE_REV_END() pInverse[id] = copy1;

__kernel void profanity_inverse_multiple(__global point * const pPoints, __global mp_number * const pInverse) {
	const size_t id = get_global_id(0) * PROFANITY_INVERSE_SIZE;

	mp_number mont_rrr = { { 0x3795f671, 0x002bb1e3, 0x00000b73, 0x1, 0, 0, 0, 0 } };
	mp_number copy1, copy2, copy3;
	mp_number b00, b01, b02, b03, b04, b05, b06, b07, b08, b09, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34, b35, b36, b37, b38, b39, b40, b41, b42, b43, b44, b45, b46, b47, b48, b49, b50, b51, b52, b53, b54, b55, b56, b57, b58, b59, b60, b61, b62, b63;
	
	INVERSE_BEGIN(b00); INVERSE(b01, b00, 1); INVERSE(b02, b01, 2); INVERSE(b03, b02, 3); INVERSE(b04, b03, 4); INVERSE(b05, b04, 5); INVERSE(b06, b05, 6); INVERSE(b07, b06, 7); INVERSE(b08, b07, 8); INVERSE(b09, b08, 9); INVERSE(b10, b09, 10); INVERSE(b11, b10, 11); INVERSE(b12, b11, 12); INVERSE(b13, b12, 13); INVERSE(b14, b13, 14); INVERSE(b15, b14, 15); INVERSE(b16, b15, 16); INVERSE(b17, b16, 17); INVERSE(b18, b17, 18); INVERSE(b19, b18, 19); INVERSE(b20, b19, 20); INVERSE(b21, b20, 21); INVERSE(b22, b21, 22); INVERSE(b23, b22, 23); INVERSE(b24, b23, 24); INVERSE(b25, b24, 25); INVERSE(b26, b25, 26); INVERSE(b27, b26, 27); INVERSE(b28, b27, 28); INVERSE(b29, b28, 29); INVERSE(b30, b29, 30); INVERSE(b31, b30, 31); INVERSE(b32, b31, 32); INVERSE(b33, b32, 33); INVERSE(b34, b33, 34); INVERSE(b35, b34, 35); INVERSE(b36, b35, 36); INVERSE(b37, b36, 37); INVERSE(b38, b37, 38); INVERSE(b39, b38, 39); INVERSE(b40, b39, 40); INVERSE(b41, b40, 41); INVERSE(b42, b41, 42); INVERSE(b43, b42, 43); INVERSE(b44, b43, 44); INVERSE(b45, b44, 45); INVERSE(b46, b45, 46); INVERSE(b47, b46, 47); INVERSE(b48, b47, 48); INVERSE(b49, b48, 49); INVERSE(b50, b49, 50); INVERSE(b51, b50, 51); INVERSE(b52, b51, 52); INVERSE(b53, b52, 53); INVERSE(b54, b53, 54); INVERSE(b55, b54, 55); INVERSE(b56, b55, 56); INVERSE(b57, b56, 57); INVERSE(b58, b57, 58); INVERSE(b59, b58, 59); INVERSE(b60, b59, 60); INVERSE(b61, b60, 61); INVERSE(b62, b61, 62); INVERSE(b63, b62, 63);

	// Middle step
	copy1 = b63;
	mp_mod_inverse(&copy1);
	mp_mul_mont(&copy1, &copy1, &mont_rrr);

	// Reverse
	INVERSE_REV(b62, 63); INVERSE_REV(b61, 62); INVERSE_REV(b60, 61); INVERSE_REV(b59, 60); INVERSE_REV(b58, 59); INVERSE_REV(b57, 58); INVERSE_REV(b56, 57); INVERSE_REV(b55, 56); INVERSE_REV(b54, 55); INVERSE_REV(b53, 54); INVERSE_REV(b52, 53); INVERSE_REV(b51, 52); INVERSE_REV(b50, 51); INVERSE_REV(b49, 50); INVERSE_REV(b48, 49); INVERSE_REV(b47, 48); INVERSE_REV(b46, 47); INVERSE_REV(b45, 46); INVERSE_REV(b44, 45); INVERSE_REV(b43, 44); INVERSE_REV(b42, 43); INVERSE_REV(b41, 42); INVERSE_REV(b40, 41); INVERSE_REV(b39, 40); INVERSE_REV(b38, 39); INVERSE_REV(b37, 38); INVERSE_REV(b36, 37); INVERSE_REV(b35, 36); INVERSE_REV(b34, 35); INVERSE_REV(b33, 34); INVERSE_REV(b32, 33); INVERSE_REV(b31, 32); INVERSE_REV(b30, 31); INVERSE_REV(b29, 30); INVERSE_REV(b28, 29); INVERSE_REV(b27, 28); INVERSE_REV(b26, 27); INVERSE_REV(b25, 26); INVERSE_REV(b24, 25); INVERSE_REV(b23, 24); INVERSE_REV(b22, 23); INVERSE_REV(b21, 22); INVERSE_REV(b20, 21); INVERSE_REV(b19, 20); INVERSE_REV(b18, 19); INVERSE_REV(b17, 18); INVERSE_REV(b16, 17); INVERSE_REV(b15, 16); INVERSE_REV(b14, 15); INVERSE_REV(b13, 14); INVERSE_REV(b12, 13); INVERSE_REV(b11, 12); INVERSE_REV(b10, 11); INVERSE_REV(b09, 10); INVERSE_REV(b08, 9); INVERSE_REV(b07, 8); INVERSE_REV(b06, 7); INVERSE_REV(b05, 6); INVERSE_REV(b04, 5); INVERSE_REV(b03, 4); INVERSE_REV(b02, 3); INVERSE_REV(b01, 2); INVERSE_REV(b00, 1); INVERSE_REV_END();
}
*/

__kernel void profanity_inverse_post(__global point * const pPoints, __global const mp_number * const pInverse) {
	const size_t id = get_global_id(0);
	
	point n = pPoints[id];
	mp_number tmp = pInverse[id];
	mp_number gx = { {0x487e2097, 0xd7362e5a, 0x29bc66db, 0x231e2953, 0x33fd129c, 0x979f48c0, 0xe9089f48, 0x9981e643} };

	// newY used as temporary variable in following two statements
	mp_mod_sub_gy( &n.y, &n.y );
	mp_mul_mont( &tmp, &tmp, &n.y );
	n.y = n.x;

	mp_mul_mont( &n.x, &tmp, &tmp );
	mp_mod_sub( &n.x, &n.x, &n.y );
	mp_mod_sub_gx(&n.x, &n.x);

	mp_mod_sub( &n.y, &gx, &n.x );
	mp_mul_mont( &n.y, &n.y, &tmp );
	mp_mod_sub_gy(&n.y, &n.y);
	
	pPoints[id] = n;
}

__kernel void profanity_end(__global point * const pPoints,	__global mp_number * const pInverse ) {
	const size_t id = get_global_id(0);
	ethhash h;
	point p = pPoints[id];

	// De-montgomerize by multiplying with one.
	mp_mul_mont_one(&p.x, &p.x);
	mp_mul_mont_one(&p.y, &p.y);

	// We can't initialize via h.q here, even if we do "h.q[i] = (ulong) 0", I have no idea why.
	for( int i = 0; i < 50; ++i ) {
		h.d[i] = 0;
	}

	h.d[0] = bswap32(p.x.d[MP_WORDS - 1]);
	h.d[1] = bswap32(p.x.d[MP_WORDS - 2]);
	h.d[2] = bswap32(p.x.d[MP_WORDS - 3]);
	h.d[3] = bswap32(p.x.d[MP_WORDS - 4]);
	h.d[4] = bswap32(p.x.d[MP_WORDS - 5]);
	h.d[5] = bswap32(p.x.d[MP_WORDS - 6]);
	h.d[6] = bswap32(p.x.d[MP_WORDS - 7]);
	h.d[7] = bswap32(p.x.d[MP_WORDS - 8]);
	h.d[8] = bswap32(p.y.d[MP_WORDS - 1]);
	h.d[9] = bswap32(p.y.d[MP_WORDS - 2]);
	h.d[10] = bswap32(p.y.d[MP_WORDS - 3]);
	h.d[11] = bswap32(p.y.d[MP_WORDS - 4]);
	h.d[12] = bswap32(p.y.d[MP_WORDS - 5]);
	h.d[13] = bswap32(p.y.d[MP_WORDS - 6]);
	h.d[14] = bswap32(p.y.d[MP_WORDS - 7]);
	h.d[15] = bswap32(p.y.d[MP_WORDS - 8]);
	h.d[16] ^= 0x01; // length 64

	sha3_keccakf(&h);

	pInverse[id].d[0] = h.d[3];
	pInverse[id].d[1] = h.d[4];
	pInverse[id].d[2] = h.d[5];
	pInverse[id].d[3] = h.d[6];
	pInverse[id].d[4] = h.d[7];
}

void profanity_result_update( const size_t id, __global const uchar * const hash, __global result * const pResult, const uchar score, const uchar scoreMax ) {
	if( score && score > scoreMax ) {
		uchar hasResult = atomic_inc(&pResult[score].found); // NOTE: If "too many" results are found it'll wrap around to 0 again and overwrite last result. Only relevant if global worksize exceeds MAX(uint).

		// Save only one result for each score, the first.
		if( hasResult == 0 ) {
			pResult[score].foundId = id;

			for( int i = 0; i < 20; ++i ) {
				pResult[score].foundHash[i] = hash[i];
			}
		}
	}
}

__kernel void profanity_transform_identity(__global mp_number * const pInverse) {
}

__kernel void profanity_transform_contract(__global mp_number * const pInverse) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;

	ethhash h;
	for( int i = 0; i < 50; ++i ) {
		h.d[i] = 0;
	}
	// set up keccak(0xd6, 0x94, address, 0x80)
	h.b[0] = 214;
	h.b[1] = 148;
	for (int i = 0; i < 20; i++) {
		h.b[i+2] = hash[i];
	}
	h.b[22] = 128;

	h.b[23] ^= 0x01; // length 23
	sha3_keccakf(&h);

	pInverse[id].d[0] = h.d[3];
	pInverse[id].d[1] = h.d[4];
	pInverse[id].d[2] = h.d[5];
	pInverse[id].d[3] = h.d[6];
	pInverse[id].d[4] = h.d[7];
}

__kernel void profanity_score_benchmark(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax ) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;
	int score = 0;

	profanity_result_update(id, hash, pResult, score, scoreMax);
}

__kernel void profanity_score_matching(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax ) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;
	int score = 0;

	for( int i = 0; i < 20; ++i ) {
		if( data1[i] > 0 && (hash[i] & data1[i]) == data2[i] ) {
			++score;
		}
	}

	profanity_result_update(id, hash, pResult, score, scoreMax);
}

__kernel void profanity_score_leading(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax ) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;
	int score = 0;

	for( int i = 0; i < 20; ++i ) {
		if( (hash[i] & 0xF0) >> 4 == data1[0] ) {
			++score;
		} else {
			break;
		}
			
		if( (hash[i] & 0x0F) == data1[0] ) {
			++score;
		} else {
			break;
		}
	}

	profanity_result_update(id, hash, pResult, score, scoreMax);
}

__kernel void profanity_score_range(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax ) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;
	int score = 0;

	for( int i = 0; i < 20; ++i ) {
		const uchar first = (hash[i] & 0xF0) >> 4;
		const uchar second = (hash[i] & 0x0F);

		if( first >= data1[0] && first <= data2[0] ) {
			++score;
		}

		if( second >= data1[0] && second <= data2[0] ) {
			++score;
		}
	}

	profanity_result_update(id, hash, pResult, score, scoreMax);
}

__kernel void profanity_score_leadingrange(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax ) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;
	int score = 0;

	for( int i = 0; i < 20; ++i ) {
		const uchar first = (hash[i] & 0xF0) >> 4;
		const uchar second = (hash[i] & 0x0F);

		if( first >= data1[0] && first <= data2[0] ) {
			++score;
		} else {
			break;
		}

		if( second >= data1[0] && second <= data2[0] ) {
			++score;
		} else {
			break;
		}
	}

	profanity_result_update(id, hash, pResult, score, scoreMax);
}

__kernel void profanity_score_mirror(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax ) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;
	int score = 0;

	for( int i = 0; i < 10; ++i ) {
		const uchar leftLeft = (hash[9-i] & 0xF0) >> 4;
		const uchar leftRight = (hash[9-i] & 0x0F); 

		const uchar rightLeft = (hash[10+i] & 0xF0) >> 4;
		const uchar rightRight = (hash[10+i] & 0x0F);

		if( leftRight != rightLeft ) {
			break;
		}

		++score;

		if( leftLeft != rightRight ) {
			break;
		}

		++score;
	}

	profanity_result_update(id, hash, pResult, score, scoreMax);
}

__kernel void profanity_score_doubles(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;
	int score = 0;

	for (int i = 0; i < 20; ++i) {
		if( (hash[i] == 0x00) || (hash[i] == 0x11) || (hash[i] == 0x22) || (hash[i] == 0x33) || (hash[i] == 0x44) || (hash[i] == 0x55) || (hash[i] == 0x66) || (hash[i] == 0x77) || (hash[i] == 0x88) || (hash[i] == 0x99) || (hash[i] == 0xAA) || (hash[i] == 0xBB) || (hash[i] == 0xCC) || (hash[i] == 0xDD) || (hash[i] == 0xEE) || (hash[i] == 0xFF) ) {
			++score;
		} else {
			break;
		}
	}

	profanity_result_update(id, hash, pResult, score, scoreMax);
}
