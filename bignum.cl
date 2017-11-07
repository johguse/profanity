/*
 * This file contains two modifications from the original:
 *   a) The file has been cut to only include the BIGNUM library
 *   b) The bignum-struct has been made into a union including a
 *      member for byte-level access.
 */

/*
 * Vanitygen, vanity bitcoin address generator
 * Copyright (C) 2011 <samr7@cs.washington.edu>
 *
 * Vanitygen is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version. 
 *
 * Vanitygen is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with Vanitygen.  If not, see <http://www.gnu.org/licenses/>.
 */

/* Byte-swapping and endianness */
#define bswap32(v)					\
	(((v) >> 24) | (((v) >> 8) & 0xff00) |		\
	 (((v) << 8) & 0xff0000) | ((v) << 24))

#if __ENDIAN_LITTLE__ != 1
#define load_le32(v) bswap32(v)
#define load_be32(v) (v)
#else
#define load_le32(v) (v)
#define load_be32(v) bswap32(v)
#endif

/*
 * Loop unrolling macros
 *
 * In most cases, preprocessor unrolling works best.
 * The exception is NVIDIA's compiler, which seems to take unreasonably
 * long to compile a loop with a larger iteration count, or a loop with
 * a body of >50 PTX instructions, with preprocessor unrolling.
 * However, it does not seem to take as long with pragma unroll, and
 * produces good output.
 */

/* Explicit loop unrolling */
#define unroll_5(a) do { a(0) a(1) a(2) a(3) a(4) } while (0)
#define unroll_8(a) do { a(0) a(1) a(2) a(3) a(4) a(5) a(6) a(7) } while (0)
#define unroll_1_7(a) do { a(1) a(2) a(3) a(4) a(5) a(6) a(7) } while (0)
#define unroll_7(a) do { a(0) a(1) a(2) a(3) a(4) a(5) a(6) } while (0)
#define unroll_7_0(a) do { a(7) a(6) a(5) a(4) a(3) a(2) a(1) a(0) } while (0)
#define unroll_7_1(a) do { a(7) a(6) a(5) a(4) a(3) a(2) a(1) } while (0)
#define unroll_16(a) do {				\
	a(0) a(1) a(2) a(3) a(4) a(5) a(6) a(7)		\
	a(8) a(9) a(10) a(11) a(12) a(13) a(14) a(15)	\
	} while (0)
#define unroll_64(a) do {				\
	a(0) a(1) a(2) a(3) a(4) a(5) a(6) a(7)		\
	a(8) a(9) a(10) a(11) a(12) a(13) a(14) a(15)	\
	a(16) a(17) a(18) a(19) a(20) a(21) a(22) a(23) \
	a(24) a(25) a(26) a(27) a(28) a(29) a(30) a(31)	\
	a(32) a(33) a(34) a(35) a(36) a(37) a(38) a(39) \
	a(40) a(41) a(42) a(43) a(44) a(45) a(46) a(47) \
	a(48) a(49) a(50) a(51) a(52) a(53) a(54) a(55) \
	a(56) a(57) a(58) a(59) a(60) a(61) a(62) a(63) \
	} while (0)

/* Conditional loop unrolling */
#if defined(DEEP_PREPROC_UNROLL)
#define iter_5(a) unroll_5(a)
#define iter_8(a) unroll_8(a)
#define iter_16(a) unroll_16(a)
#define iter_64(a) unroll_64(a)
#else
#define iter_5(a) do {int _i; for (_i = 0; _i < 5; _i++) { a(_i) }} while (0)
#define iter_8(a) do {int _i; for (_i = 0; _i < 8; _i++) { a(_i) }} while (0)
#define iter_16(a) do {int _i; for (_i = 0; _i < 16; _i++) { a(_i) }} while (0)
#define iter_64(a) do {int _i; for (_i = 0; _i < 64; _i++) { a(_i) }} while (0)
#endif

/*
 * BIGNUM mini-library
 * This module deals with fixed-size 256-bit bignums.
 * Where modular arithmetic is performed, the SECP256k1 prime
 * modulus (below) is assumed.
 *
 * Methods include:
 * - bn_is_zero/bn_is_one/bn_is_odd/bn_is_even/bn_is_bit_set
 * - bn_rshift[1]/bn_lshift[1]
 * - bn_neg
 * - bn_uadd/bn_uadd_p
 * - bn_usub/bn_usub_p
 */

typedef uint bn_word;
#define BN_NBITS 256
#define BN_WSHIFT 5
#define BN_WBITS (1 << BN_WSHIFT)
#define BN_NWORDS ((BN_NBITS/8) / sizeof(bn_word))
#define BN_WORDMAX 0xffffffff

#define MODULUS_BYTES \
	0xfffffc2f, 0xfffffffe, 0xffffffff, 0xffffffff, \
	0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff

typedef union {
	bn_word d[BN_NWORDS];
	uchar b[BN_NWORDS * 4];
} bignum;

__constant bn_word modulus[] = { MODULUS_BYTES };
__constant bignum bn_zero = { 0 };

__constant bn_word mont_rr[BN_NWORDS] = { 0xe90a1, 0x7a2, 0x1, 0, };
__constant bn_word mont_n0[2] = { 0xd2253531, 0xd838091d };


#define bn_is_odd(bn)		(bn.d[0] & 1)
#define bn_is_even(bn) 		(!bn_is_odd(bn))
#define bn_is_zero(bn) 		(!bn.d[0] && !bn.d[1] && !bn.d[2] && \
				 !bn.d[3] && !bn.d[4] && !bn.d[5] && \
				 !bn.d[6] && !bn.d[7])
#define bn_is_one(bn) 		((bn.d[0] == 1) && !bn.d[1] && !bn.d[2] && \
				 !bn.d[3] && !bn.d[4] && !bn.d[5] && \
				 !bn.d[6] && !bn.d[7])
#define bn_is_bit_set(bn, n) \
	((((bn_word*)&bn)[n >> BN_WSHIFT]) & (1 << (n & (BN_WBITS-1))))

#define bn_unroll(e) unroll_8(e)
#define bn_unroll_sf(e)	unroll_1_7(e)
#define bn_unroll_sl(e)	unroll_7(e)
#define bn_unroll_reverse(e) unroll_7_0(e)
#define bn_unroll_reverse_sl(e) unroll_7_1(e)

#define bn_unroll_arg(e, arg)				\
	e(arg, 0) e(arg, 1) e(arg, 2) e(arg, 3)	\
	e(arg, 4) e(arg, 5) e(arg, 6) e(arg, 7)
#define bn_unroll_arg_sf(e, arg)			\
	e(arg, 1) e(arg, 2) e(arg, 3)		\
	e(arg, 4) e(arg, 5) e(arg, 6) e(arg, 7)

#define bn_iter(e) iter_8(e)


/*
 * Bitwise shift
 */

void
bn_lshift1(bignum *bn)
{
#define bn_lshift1_inner1(i)						\
		bn->d[i] = (bn->d[i] << 1) | (bn->d[i-1] >> 31);
	bn_unroll_reverse_sl(bn_lshift1_inner1);
	bn->d[0] <<= 1;
}

void
bn_rshift(bignum *bn, int shift)
{
	int wd, iws, iwr;
	bn_word ihw, ilw;
	iws = (shift & (BN_WBITS-1));
	iwr = BN_WBITS - iws;
	wd = (shift >> BN_WSHIFT);
	ihw = (wd < BN_WBITS) ? bn->d[wd] : 0;

#define bn_rshift_inner1(i)				\
		wd++;					\
		ilw = ihw;				\
		ihw = (wd < BN_WBITS) ? bn->d[wd] : 0;	\
		bn->d[i] = (ilw >> iws) | (ihw << iwr);
	bn_unroll_sl(bn_rshift_inner1);
	bn->d[BN_NWORDS-1] = (ihw >> iws);
}

void
bn_rshift1(bignum *bn)
{
#define bn_rshift1_inner1(i)						\
		bn->d[i] = (bn->d[i+1] << 31) | (bn->d[i] >> 1);
	bn_unroll_sl(bn_rshift1_inner1);
	bn->d[BN_NWORDS-1] >>= 1;
}

void
bn_rshift1_2(bignum *bna, bignum *bnb)
{
#define bn_rshift1_2_inner1(i)						\
		bna->d[i] = (bna->d[i+1] << 31) | (bna->d[i] >> 1);	\
		bnb->d[i] = (bnb->d[i+1] << 31) | (bnb->d[i] >> 1);
	bn_unroll_sl(bn_rshift1_2_inner1);
	bna->d[BN_NWORDS-1] >>= 1;
	bnb->d[BN_NWORDS-1] >>= 1;
}


/*
 * Unsigned comparison
 */

int
bn_ucmp_ge(bignum *a, bignum *b)
{
	int l = 0, g = 0;

#define bn_ucmp_ge_inner1(i)				\
		if (a->d[i] < b->d[i]) l |= (1 << i);	\
		if (a->d[i] > b->d[i]) g |= (1 << i);
	bn_unroll_reverse(bn_ucmp_ge_inner1);
	return (l > g) ? 0 : 1;
}

int
bn_ucmp_ge_c(bignum *a, __constant bn_word *b)
{
	int l = 0, g = 0;

#define bn_ucmp_ge_c_inner1(i)				\
		if (a->d[i] < b[i]) l |= (1 << i);	\
		if (a->d[i] > b[i]) g |= (1 << i);
	bn_unroll_reverse(bn_ucmp_ge_c_inner1);
	return (l > g) ? 0 : 1;
}

/*
 * Negate
 */

void
bn_neg(bignum *n)
{
	int c = 1;

#define bn_neg_inner1(i)				\
		c = (n->d[i] = (~n->d[i]) + c) ? 0 : c;
	bn_unroll(bn_neg_inner1);
}

/*
 * Add/subtract
 */

#define bn_add_word(r, a, b, t, c) do {		\
		t = a + b;			\
		c = (t < a) ? 1 : 0;		\
		r = t;				\
	} while (0)

#define bn_addc_word(r, a, b, t, c) do {			\
		t = a + b + c;					\
		c = (t < a) ? 1 : ((c & (t == a)) ? 1 : 0);	\
		r = t;						\
	} while (0)

bn_word
bn_uadd_words_seq(bn_word *r, bn_word *a, bn_word *b)
{
	bn_word t, c = 0;

#define bn_uadd_words_seq_inner1(i)			\
		bn_addc_word(r[i], a[i], b[i], t, c);
	bn_add_word(r[0], a[0], b[0], t, c);
	bn_unroll_sf(bn_uadd_words_seq_inner1);
	return c;
}

bn_word
bn_uadd_words_c_seq(bn_word *r, bn_word *a, __constant bn_word *b)
{
	bn_word t, c = 0;

	bn_add_word(r[0], a[0], b[0], t, c);
	bn_unroll_sf(bn_uadd_words_seq_inner1);
	return c;
}

#define bn_sub_word(r, a, b, t, c) do {		\
		t = a - b;			\
		c = (a < b) ? 1 : 0;		\
		r = t;				\
	} while (0)

#define bn_subb_word(r, a, b, t, c) do {	\
		t = a - (b + c);		\
		c = (!(a) && c) ? 1 : 0;	\
		c |= (a < b) ? 1 : 0;		\
		r = t;				\
	} while (0)

bn_word
bn_usub_words_seq(bn_word *r, bn_word *a, bn_word *b)
{
	bn_word t, c = 0;

#define bn_usub_words_seq_inner1(i)			\
		bn_subb_word(r[i], a[i], b[i], t, c);

	bn_sub_word(r[0], a[0], b[0], t, c);
	bn_unroll_sf(bn_usub_words_seq_inner1);
	return c;
}

bn_word
bn_usub_words_c_seq(bn_word *r, bn_word *a, __constant bn_word *b)
{
	bn_word t, c = 0;

	bn_sub_word(r[0], a[0], b[0], t, c);
	bn_unroll_sf(bn_usub_words_seq_inner1);
	return c;
}

/*
 * Add/subtract better suited for AMD's VLIW architecture
 */
bn_word
bn_uadd_words_vliw(bn_word *r, bn_word *a, bn_word *b)
{
	bignum x;
	bn_word c = 0, cp = 0;

#define bn_uadd_words_vliw_inner1(i)		\
		x.d[i] = a[i] + b[i];

#define bn_uadd_words_vliw_inner2(i)			\
		c |= (a[i] > x.d[i]) ? (1 << i) : 0;	\
		cp |= (!~x.d[i]) ? (1 << i) : 0;

#define bn_uadd_words_vliw_inner3(i)		\
		r[i] = x.d[i] + ((c >> i) & 1);

	bn_unroll(bn_uadd_words_vliw_inner1);
	bn_unroll(bn_uadd_words_vliw_inner2);
	c = ((cp + (c << 1)) ^ cp);
	r[0] = x.d[0];
	bn_unroll_sf(bn_uadd_words_vliw_inner3);
	return c >> BN_NWORDS;
}

bn_word
bn_uadd_words_c_vliw(bn_word *r, bn_word *a, __constant bn_word *b)
{
	bignum x;
	bn_word c = 0, cp = 0;

	bn_unroll(bn_uadd_words_vliw_inner1);
	bn_unroll(bn_uadd_words_vliw_inner2);
	c = ((cp + (c << 1)) ^ cp);
	r[0] = x.d[0];
	bn_unroll_sf(bn_uadd_words_vliw_inner3);
	return c >> BN_NWORDS;
}

bn_word
bn_usub_words_vliw(bn_word *r, bn_word *a, bn_word *b)
{
	bignum x;
	bn_word c = 0, cp = 0;

#define bn_usub_words_vliw_inner1(i)		\
		x.d[i] = a[i] - b[i];

#define bn_usub_words_vliw_inner2(i)			\
		c |= (a[i] < b[i]) ? (1 << i) : 0;	\
		cp |= (!x.d[i]) ? (1 << i) : 0;

#define bn_usub_words_vliw_inner3(i)		\
		r[i] = x.d[i] - ((c >> i) & 1);

	bn_unroll(bn_usub_words_vliw_inner1);
	bn_unroll(bn_usub_words_vliw_inner2);
	c = ((cp + (c << 1)) ^ cp);
	r[0] = x.d[0];
	bn_unroll_sf(bn_usub_words_vliw_inner3);
	return c >> BN_NWORDS;
}

bn_word
bn_usub_words_c_vliw(bn_word *r, bn_word *a, __constant bn_word *b)
{
	bignum x;
	bn_word c = 0, cp = 0;

	bn_unroll(bn_usub_words_vliw_inner1);
	bn_unroll(bn_usub_words_vliw_inner2);
	c = ((cp + (c << 1)) ^ cp);
	r[0] = x.d[0];
	bn_unroll_sf(bn_usub_words_vliw_inner3);
	return c >> BN_NWORDS;
}


#if defined(DEEP_VLIW)
#define bn_uadd_words bn_uadd_words_vliw
#define bn_uadd_words_c bn_uadd_words_c_vliw
#define bn_usub_words bn_usub_words_vliw
#define bn_usub_words_c bn_usub_words_c_vliw
#else
#define bn_uadd_words bn_uadd_words_seq
#define bn_uadd_words_c bn_uadd_words_c_seq
#define bn_usub_words bn_usub_words_seq
#define bn_usub_words_c bn_usub_words_c_seq
#endif

#define bn_uadd(r, a, b) bn_uadd_words((r)->d, (a)->d, (b)->d)
#define bn_uadd_c(r, a, b) bn_uadd_words_c((r)->d, (a)->d, b)
#define bn_usub(r, a, b) bn_usub_words((r)->d, (a)->d, (b)->d)
#define bn_usub_c(r, a, b) bn_usub_words_c((r)->d, (a)->d, b)

/*
 * Modular add/sub
 */

void
bn_mod_add(bignum *r, bignum *a, bignum *b)
{
	if (bn_uadd(r, a, b) ||
	    (bn_ucmp_ge_c(r, modulus)))
		bn_usub_c(r, r, modulus);
}

void
bn_mod_sub(bignum *r, bignum *a, bignum *b)
{
	if (bn_usub(r, a, b))
		bn_uadd_c(r, r, modulus);
}

void
bn_mod_lshift1(bignum *bn)
{
	bn_word c = (bn->d[BN_NWORDS-1] & 0x80000000);
	bn_lshift1(bn);
	if (c || (bn_ucmp_ge_c(bn, modulus)))
		bn_usub_c(bn, bn, modulus);
}

/*
 * Montgomery multiplication
 *
 * This includes normal multiplication of two "Montgomeryized"
 * bignums, and bn_from_mont for de-Montgomeryizing a bignum.
 */

#define bn_mul_word(r, a, w, c, p, s) do { \
		r = (a * w) + c;	   \
		p = mul_hi(a, w);	   \
		c = (r < c) ? p + 1 : p;   \
	} while (0)

#define bn_mul_add_word(r, a, w, c, p, s) do {	\
		s = r + c;			\
		p = mul_hi(a, w);		\
		r = (a * w) + s;		\
		c = (s < c) ? p + 1 : p;	\
		if (r < s) c++;			\
	} while (0)
void
bn_mul_mont(bignum *r, bignum *a, bignum *b)
{
	bignum t;
	bn_word tea, teb, c, p, s, m;

#if !defined(VERY_EXPENSIVE_BRANCHES)
	int q;
#endif

	c = 0;
#define bn_mul_mont_inner1(j)					\
		bn_mul_word(t.d[j], a->d[j], b->d[0], c, p, s);
	bn_unroll(bn_mul_mont_inner1);
	tea = c;
	teb = 0;

	c = 0;
	m = t.d[0] * mont_n0[0];
	bn_mul_add_word(t.d[0], modulus[0], m, c, p, s);
#define bn_mul_mont_inner2(j)						\
		bn_mul_add_word(t.d[j], modulus[j], m, c, p, s);	\
		t.d[j-1] = t.d[j];
	bn_unroll_sf(bn_mul_mont_inner2);
	t.d[BN_NWORDS-1] = tea + c;
	tea = teb + ((t.d[BN_NWORDS-1] < c) ? 1 : 0);

#define bn_mul_mont_inner3_1(i, j)					\
		bn_mul_add_word(t.d[j], a->d[j], b->d[i], c, p, s);
#define bn_mul_mont_inner3_2(i, j)					\
		bn_mul_add_word(t.d[j], modulus[j], m, c, p, s);	\
		t.d[j-1] = t.d[j];
#define bn_mul_mont_inner3(i)				 \
	c = 0;						 \
	bn_unroll_arg(bn_mul_mont_inner3_1, i);		 \
	tea += c;					 \
	teb = ((tea < c) ? 1 : 0);			 \
	c = 0;						 \
	m = t.d[0] * mont_n0[0];			 \
	bn_mul_add_word(t.d[0], modulus[0], m, c, p, s); \
	bn_unroll_arg_sf(bn_mul_mont_inner3_2, i);	 \
	t.d[BN_NWORDS-1] = tea + c;			 \
	tea = teb + ((t.d[BN_NWORDS-1] < c) ? 1 : 0);

	/*
	 * The outer loop here is quite long, and we won't unroll it
	 * unless VERY_EXPENSIVE_BRANCHES is set.
	 */
#if defined(VERY_EXPENSIVE_BRANCHES)
	bn_unroll_sf(bn_mul_mont_inner3);
	c = tea | !bn_usub_c(r, &t, modulus);
	if (!c)
		*r = t;

#else
	for (q = 1; q < BN_NWORDS; q++) {
		bn_mul_mont_inner3(q);
	}
	c = tea || (t.d[BN_NWORDS-1] >= modulus[BN_NWORDS-1]);
	if (c) {
		c = tea | !bn_usub_c(r, &t, modulus);
		if (c)
			return;
	}
	*r = t;
#endif
}

void
bn_from_mont(bignum *rb, bignum *b)
{
#define WORKSIZE ((2*BN_NWORDS) + 1)
	bn_word r[WORKSIZE];
	bn_word m, c, p, s;
#if defined(PRAGMA_UNROLL)
	int i;
#endif

	/* Copy the input to the working area */
	/* Zero the upper words */
#define bn_from_mont_inner1(i)			\
	r[i] = b->d[i];
#define bn_from_mont_inner2(i)			\
	r[BN_NWORDS+i] = 0;

	bn_unroll(bn_from_mont_inner1);
	bn_unroll(bn_from_mont_inner2);
	r[WORKSIZE-1] = 0;

	/* Multiply (long) by modulus */
#define bn_from_mont_inner3_1(i, j) \
	bn_mul_add_word(r[i+j], modulus[j], m, c, p, s);

#if !defined(VERY_EXPENSIVE_BRANCHES)
#define bn_from_mont_inner3_2(i)		\
	if (r[BN_NWORDS + i] < c)		\
		r[BN_NWORDS + i + 1] += 1;
#else
#define bn_from_mont_inner3_2(i)				\
	r[BN_NWORDS + i + 1] += (r[BN_NWORDS + i] < c) ? 1 : 0;
#endif

#define bn_from_mont_inner3(i)			 \
	m = r[i] * mont_n0[0];			 \
	c = 0;					 \
	bn_unroll_arg(bn_from_mont_inner3_1, i); \
	r[BN_NWORDS + i] += c;			 \
	bn_from_mont_inner3_2(i)

	/*
	 * The outer loop here is not very long, so we will unroll
	 * it by default.  However, it's just complicated enough to
	 * cause NVIDIA's compiler to take unreasonably long to compile
	 * it, unless we use pragma unroll.
	 */
#if !defined(PRAGMA_UNROLL)
	bn_iter(bn_from_mont_inner3);
#else
#pragma unroll 8
	for (i = 0; i < BN_NWORDS; i++) { bn_from_mont_inner3(i) }
#endif

	/*
	 * Make sure the result is less than the modulus.
	 * Subtracting is not much more expensive than compare, so
	 * subtract always and assign based on the carry out value.
	 */
	c = bn_usub_words_c(rb->d, &r[BN_NWORDS], modulus);
	if (c) {
#define bn_from_mont_inner4(i)				\
			rb->d[i] = r[BN_NWORDS + i];
		bn_unroll(bn_from_mont_inner4);
	}
}

/*
 * Modular inversion
 */

void
bn_mod_inverse(bignum *r, bignum *n)
{
	bignum a, b, x, y;
	int shift;
	bn_word xc, yc;
	for (shift = 0; shift < BN_NWORDS; shift++) {
		a.d[shift] = modulus[shift];
		x.d[shift] = 0;
		y.d[shift] = 0;
	}
	b = *n;
	x.d[0] = 1;
	xc = 0;
	yc = 0;
	while (!bn_is_zero(b)) {
		shift = 0;
		while (!bn_is_odd(b)) {
			if (bn_is_odd(x))
				xc += bn_uadd_c(&x, &x, modulus);
			bn_rshift1_2(&x, &b);
			x.d[7] |= (xc << 31);
			xc >>= 1;
		}

		while (!bn_is_odd(a)) {
			if (bn_is_odd(y))
				yc += bn_uadd_c(&y, &y, modulus);
			bn_rshift1_2(&y, &a);
			y.d[7] |= (yc << 31);
			yc >>= 1;
		}

		if (bn_ucmp_ge(&b, &a)) {
			xc += yc + bn_uadd(&x, &x, &y);
			bn_usub(&b, &b, &a);
		} else {
			yc += xc + bn_uadd(&y, &y, &x);
			bn_usub(&a, &a, &b);
		}
	}

	if (!bn_is_one(a)) {
		/* no modular inverse */
		*r = bn_zero;
	} else {
		/* Compute y % m as cheaply as possible */
		while (yc < 0x80000000)
			yc -= bn_usub_c(&y, &y, modulus);
		bn_neg(&y);
		*r = y;
	}
}
