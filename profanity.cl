/* ------------------------------------------------------------------------ */
/* Elliptic point and addition (with caveats).                              */
/* ------------------------------------------------------------------------ */
typedef struct {
    bignum x;
    bignum y;
} point;

// Does not handle:
//   * Points sharing X coordinate
void point_add(point * const p, point * const o) {
	bignum monter = { 0xe90a1, 0x7a2, 0x1, 0 };
	bignum deltaY;
	bignum tmp;
	bignum newX;
	bignum newY;

	bn_mod_sub( &deltaY, &o->y, &p->y );

	bn_mod_sub( &tmp, &o->x, &p->x );
	bn_mod_inverse( &tmp, &tmp );
	bn_mul_mont( &tmp, &tmp, &deltaY );
	bn_mul_mont( &tmp, &tmp, &monter );

	bn_mul_mont( &newX, &tmp, &tmp );
	bn_mul_mont( &newX, &newX, &monter );
	bn_mod_sub( &newX, &newX, &p->x );
	bn_mod_sub( &newX, &newX, &o->x );

	bn_mod_sub( &newY, &p->x, &newX );
	bn_mul_mont( &newY, &newY, &tmp );
	bn_mul_mont( &newY, &newY, &monter );
	bn_mod_sub( &newY, &newY, &p->y );

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
		pResult->found = 0;
	}
}

__kernel void profanity_inverse_pre(__global const point * const precomp, __global const point * const pPoints1, __global point * const pPoints2, __global bignum * const pInverse, __global uchar * pPass ) {
	const size_t id = get_global_id(0);
	
	point s = pPoints1[id / 255];
	point o = precomp[8 * 255 * 3 + (*pPass) * 255 + id % 255];
	bignum deltaX;
	
	bn_mod_sub( &deltaX, &o.x, &s.x);
	pInverse[id] = deltaX;
	pPoints2[id / 255] = s; // Multiple overwrites
}


__kernel void profanity_inverse(__global bignum * const pInverse, __global uchar * pPass ) {
	const size_t id = get_global_id(0) * 3;
	
	bignum monter = { 0xe90a1, 0x7a2, 0x1, 0 };
	bignum a = pInverse[id+0];
	bignum b = pInverse[id+1];
	bignum c = pInverse[id+2];

	bignum ab;
	bignum bc;
	bignum ca;

	bignum i;

	bn_mul_mont(&ab, &a, &b);
	bn_mul_mont(&bc, &b, &c);
	bn_mul_mont(&ca, &c, &a);

	bn_mul_mont(&i, &ab, &c);
	bn_mul_mont(&i, &i, &monter);
	bn_mod_inverse( &i, &i );

	bn_mul_mont(&ab, &ab, &i);
	bn_mul_mont(&bc, &bc, &i);
	bn_mul_mont(&ca, &ca, &i);

	bn_mul_mont(&ab, &ab, &monter);
	bn_mul_mont(&bc, &bc, &monter);
	bn_mul_mont(&ca, &ca, &monter);

	pInverse[id+0] = bc;
	pInverse[id+1] = ca;
	pInverse[id+2] = ab;
	
	// We increase the pass counter here where it's not used. (*pPass - 1) used in profanity_inverse_post
	if( id == 0 ) {
		*pPass += 1;
	}
}

__kernel void profanity_inverse_post(__global const point * const precomp, __global point * const pPoints1, __global point * const pPoints2, __global const bignum * const pInverse, __global uchar * pPass ) {
	const size_t id = get_global_id(0);
	
	point s = pPoints2[id / 255];
	point o = precomp[8 * 255 * 3 + (*pPass - 1) * 255 + id % 255];
	
	bignum monter = { 0xe90a1, 0x7a2, 0x1, 0 };
	
	bignum deltaY;
	bignum tmp = pInverse[id];
	bignum newX;
	bignum newY;

	bn_mod_sub( &deltaY, &o.y, &s.y );

	bn_mul_mont( &tmp, &tmp, &deltaY );
	bn_mul_mont( &tmp, &tmp, &monter );

	bn_mul_mont( &newX, &tmp, &tmp );
	bn_mul_mont( &newX, &newX, &monter );
	bn_mod_sub( &newX, &newX, &s.x );
	bn_mod_sub( &newX, &newX, &o.x );

	bn_mod_sub( &newY, &s.x, &newX );
	bn_mul_mont( &newY, &newY, &tmp );
	bn_mul_mont( &newY, &newY, &monter );
	bn_mod_sub( &newY, &newY, &s.y );
	
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
	__global const point * const pSelf = &pPoints1[id];
	
	ethhash h;
	
	for( uchar i = 0; i < 200; ++i ) {
		if( i < 32 ) {
			h.b[i] = pSelf->x.b[31-i];
		} else if( i < 64 ) {
			h.b[i] = pSelf->y.b[63-i];
		} else {
			h.b[i] = 0;
		}
	}
		
	sha3_keccakf(&h);
	uchar * hash = h.b + 12;
	
	uchar score = 0;
	if( mode == 0 ) {
		for( uchar i = 0; i < 20; ++i ) {
			if( (hash[i] & 0xF0) == 0 ) {
				++score;
			}
			
			if( (hash[i] & 0x0F) == 0 ) {
				++score;
			}
		}
	} else if( mode == 1 ) {
		for( uchar i = 0; i < 20; ++i ) {
			if( data1[i] > 0 && (hash[i] & data1[i]) == data2[i] ) {
				++score;
			}
		}
	} else if( mode == 2 ) {
		for( uchar i = 0; i < 20; ++i ) {
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
		for( uchar i = 0; i < 20; ++i ) {
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
		for( uchar i = 0; i < 20; ++i ) {
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
		uchar hasResult = atomic_inc(&pResult->found); // NOTE: If "too many" results are found it'll wrap around to 0 again and overwrite last result. Only relevant if global worksize exceeds MAX(uint).

		// Save only one result, the first.
		if( hasResult == 0 ) {
			pResult->foundId = id;
			pResult->foundScore = score;

			for( uchar i = 0; i < 20; ++i ) {
				pResult->foundHash[i] = hash[i];
			}
		}
	}	
}