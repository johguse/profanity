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
	uchar foundScore;
	uchar foundHash[20];
} result;

void profanity_begin_seed(__global const point * const precomp, point * const p, bool * const pIsFirst, const uchar byteCount, const size_t precompOffset, const ulong seed) {
	point o;

	for( uchar i = 0; i < byteCount; ++i ) {
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

__kernel void profanity_begin(__global const point * const precomp, __global point * const pPoints1, __global uchar * const pPass, __global result * const pResult, const ulong4 seed) {
	const size_t id = get_global_id(0);

	if( id == 0 ) {
		point p;
		point o;
		bool bIsFirst = true;

		profanity_begin_seed(precomp, &p, &bIsFirst, 8, 8 * 255 * 0, seed.x);
		profanity_begin_seed(precomp, &p, &bIsFirst, 8, 8 * 255 * 1, seed.y);
		profanity_begin_seed(precomp, &p, &bIsFirst, 8, 8 * 255 * 2, seed.z);
		profanity_begin_seed(precomp, &p, &bIsFirst, 8 - PROFANITY_PASSES, 8 * 255 * 3, seed.w);

		pPoints1[0] = p;
		*pPass = 8 - PROFANITY_PASSES;
		for( uchar i = 0; i < 40; ++i ) {
			pResult[i].found = 0;
		}
	}
}

__kernel void profanity_inverse_pre(__global const point * const precomp, __global const point * const pPoints1, __global point * const pPoints2, __global mp_number * const pInverse, __global uchar * pPass ) {
	const size_t id = get_global_id(0);
	
	point s = pPoints1[id / 255];
	point o = precomp[8 * 255 * 3 + (*pPass) * 255 + id % 255];
	mp_number deltaX;
	
	mp_mod_sub( &deltaX, &o.x, &s.x);
	pInverse[id] = deltaX;
	pPoints2[id / 255] = s; // Multiple overwrites
}

__kernel void profanity_inverse_multiple(__global mp_number * const pInverse, __global uchar * pPass ) {
	const size_t id = get_global_id(0) * 255;
	
	mp_number inv;
	mp_number copy; // Optimize this later
	mp_number buffer[255];
	mp_number mont_rrr = { { 0x3795f671, 0x002bb1e3, 0x00000b73, 0x1, 0, 0, 0, 0 } };
	
	buffer[0] = pInverse[id];
	for( uchar i = 1; i < 255; ++i ) {
		copy = pInverse[id + i];
		mp_mul_mont( &buffer[i], &buffer[i-1], &copy );
	}

	// mp_mod_inverse(aR) -> (aR)^-1 
	// mp_mul_mont(x,y) -> x * y * R^-1
	// mp_mul_mont( (aR)^-1, R^3 ) -> (aR)^-1 * R^3 * R^-1 = a^-1 * R
	// Also: Compiler really fucks the below up unless we use a temporary variable. Why?!
	inv = buffer[255-1];
	mp_mod_inverse( &inv );
	mp_mul_mont(&inv, &inv, &mont_rrr);

	for( uchar i = 255 - 1; i > 0; --i ) {
		copy = pInverse[id+i];
		mp_mul_mont( &copy, &copy, &inv);

		mp_mul_mont( &buffer[i], &buffer[i-1], &inv);
		pInverse[id+i] = buffer[i];
		inv = copy;
	}

	pInverse[id] = inv;

	// We increase the pass counter here where it's not used. (*pPass - 1) used in profanity_inverse_post
	if( id == 0 ) {
		*pPass += 1;
	}
}

__kernel void profanity_inverse_post(__global const point * const precomp, __global point * const pPoints1, __global point * const pPoints2, __global const mp_number * const pInverse, __global uchar * pPass ) {
	const size_t id = get_global_id(0);
	
	point s = pPoints2[id / 255];
	point o = precomp[8 * 255 * 3 + (*pPass - 1) * 255 + id % 255];
	
	mp_number tmp = pInverse[id];
	mp_number newX;
	mp_number newY;

	// newX used as temporary variable in following two statements
	mp_mod_sub( &newX, &o.y, &s.y );
	mp_mul_mont( &tmp, &tmp, &newX );

	mp_mul_mont( &newX, &tmp, &tmp );
	mp_mod_sub( &newX, &newX, &s.x );
	mp_mod_sub( &newX, &newX, &o.x );

	mp_mod_sub( &newY, &s.x, &newX );
	mp_mul_mont( &newY, &newY, &tmp );
	mp_mod_sub( &newY, &newY, &s.y );
	
	pPoints1[id].x = newX;
	pPoints1[id].y = newY;
}

__kernel void profanity_end(
	__global point * const pPoints1,
	__global result * const pResult,
	__constant const uchar * const data1,
	__constant const uchar * const data2,
	const uchar scoreMax,
	const uchar mode )
{
	const size_t id = get_global_id(0);
	ethhash h = { { 0 } }; // This doesn't work for some reason, we zero-initialize below.
	point self = pPoints1[id];
	uchar i;

	// De-montgomerize by multiplying with one.
	mp_mul_mont_one(&self.x, &self.x);
	mp_mul_mont_one(&self.y, &self.y);

	// We can't initialize via h.q here, even if we do "h.q[i] = (ulong) 0", I have no idea why.
	for( i = 0; i < 50; ++i ) {
		h.d[i] = 0;
	}

	for( i = 0; i < MP_WORDS; ++i ) {
		h.d[i] = bswap32( self.x.d[MP_WORDS - 1 - i] );
		h.d[i+8] = bswap32( self.y.d[MP_WORDS - 1 - i] );
	}
	
	sha3_keccakf(&h);
	const uchar * const hash = h.b + 12;
	uchar score = 0;
	
	if( mode == 0 ) {
		for( i = 0; i < 20; ++i ) {
			if( (hash[i] & 0xF0) == 0 ) {
				++score;
			}
			
			if( (hash[i] & 0x0F) == 0 ) {
				++score;
			}
		}
	} else if( mode == 1 ) {
		for( i = 0; i < 20; ++i ) {
			if( data1[i] > 0 && (hash[i] & data1[i]) == data2[i] ) {
				++score;
			}
		}
	} else if( mode == 2 ) {
		for( i = 0; i < 20; ++i ) {
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
	} else if( mode == 3 ) {
		for( i = 0; i < 20; ++i ) {
			const uchar first = (hash[i] & 0xF0) >> 4;
			const uchar second = (hash[i] & 0x0F);

			if( first >= data1[0] && first <= data2[0] ) {
				++score;
			}

			if( second >= data1[0] && second <= data2[0] ) {
				++score;
			}
		}	
	} else if( mode == 4 ) {
		for( i = 0; i < 20; ++i ) {
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
	}

	if( score && score > scoreMax ) {
		uchar hasResult = atomic_inc(&pResult[score].found); // NOTE: If "too many" results are found it'll wrap around to 0 again and overwrite last result. Only relevant if global worksize exceeds MAX(uint).

		// Save only one result for each score, the first.
		if( hasResult == 0 ) {
			pResult[score].foundId = id;
			pResult[score].foundScore = score;

			for( i = 0; i < 20; ++i ) {
				pResult[score].foundHash[i] = hash[i];
			}
		}
	}	
}