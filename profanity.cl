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

	for( uchar i = 0; i < 40; ++i ) {
		pResult[i].found = 0;
	}
}

__kernel void profanity_inverse_multiple(__global point * const pPoints, __global mp_number * const pInverse) {
	const size_t id = get_global_id(0) * PROFANITY_INVERSE_SIZE;
	
	mp_number inv;
	mp_number copy; // Optimize this later
	mp_number buffer[PROFANITY_INVERSE_SIZE];
	mp_number mont_rrr = { { 0x3795f671, 0x002bb1e3, 0x00000b73, 0x1, 0, 0, 0, 0 } };

	copy = pPoints[id].x;
	mp_mod_sub_gx( &copy, &copy );
	buffer[0] = copy;
	pInverse[0] = copy;

	for( uint i = 1; i < PROFANITY_INVERSE_SIZE; ++i ) {
		copy = pPoints[id+i].x;
		mp_mod_sub_gx( &copy, &copy );
		pInverse[id+i] = copy;
		mp_mul_mont( &buffer[i], &buffer[i-1], &copy );
	}

	// mp_mod_inverse(aR) -> (aR)^-1 
	// mp_mul_mont(x,y) -> x * y * R^-1
	// mp_mul_mont( (aR)^-1, R^3 ) -> (aR)^-1 * R^3 * R^-1 = a^-1 * R
	// Also: Compiler really fucks the below up unless we use a temporary variable. Why?!
	inv = buffer[PROFANITY_INVERSE_SIZE-1];
	mp_mod_inverse( &inv );
	mp_mul_mont(&inv, &inv, &mont_rrr); 

	for( uint i = PROFANITY_INVERSE_SIZE - 1; i > 0; --i ) {
		copy = pInverse[id+i];
		mp_mul_mont( &copy, &copy, &inv);

		mp_mul_mont( &buffer[i], &buffer[i-1], &inv);
		pInverse[id+i] = buffer[i];
		inv = copy;
	}

	pInverse[id] = inv;
}

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

__kernel void profanity_end(
	__global point * const pPoints,
	__global mp_number * const pInverse ) {
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

	for( int i = 0; i < MP_WORDS; ++i ) {
		h.d[i] = bswap32( p.x.d[MP_WORDS - 1 - i] );
		h.d[i+8] = bswap32( p.y.d[MP_WORDS - 1 - i] );
	}
	
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