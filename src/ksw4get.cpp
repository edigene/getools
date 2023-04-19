#include "ksw4get.h"

/**
 * Initialize the query data structure
 *
 * @param size   Number of bytes used to store a score; valid valures are 1 or 2
 * @param qlen   Length of the query sequence
 * @param query  Query sequence
 * @param m      Size of the alphabet
 * @param mat    Scoring matrix in a one-dimension array
 *
 * @return       Query data structure
 */
kswq_t *kswge_qinit(int size, int qlen, const uint8_t *query, int m, const int8_t *mat)
{
	kswq_t *q;
	int slen, a, tmp, p;

	size = size > 1? 2 : 1;
	p = 8 * (3 - size); // # values per __m128i
	slen = (qlen + p - 1) / p; // segmented length
	q = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (m + 4)); // a single block of memory
	q->qp = (__m128i*)(((size_t)q + sizeof(kswq_t) + 15) >> 4 << 4); // align memory
	q->H0 = q->qp + slen * m;
	q->H1 = q->H0 + slen;
	q->E  = q->H1 + slen;
	q->Hmax = q->E + slen;
	q->slen = slen; q->qlen = qlen; q->size = size;
	// compute shift
	tmp = m * m;
	for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
		if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
		if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
	}
	q->max = q->mdiff;
	q->shift = 256 - q->shift; // NB: q->shift is uint8_t
	q->mdiff += q->shift; // this is the difference between the min and max scores
	// An example: p=8, qlen=19, slen=3 and segmentation:
	//  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
	if (size == 1) {
		int8_t *t = (int8_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]) + q->shift;
		}
	} else {
		int16_t *t = (int16_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]);
		}
	}
	return q;
}

kswr_t kswge_u8(kswq_t *q, int tlen, const uint8_t *target, int _o_del, int _e_del, int _o_ins, int _e_ins, int xtra) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i zero, oe_del, e_del, oe_ins, e_ins, shift, *H0, *H1, *E, *Hmax;
	kswr_t r;

#define __max_16(ret, xx) do { \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
    	(ret) = _mm_extract_epi16((xx), 0) & 0x00ff; \
	} while (0)

	// initialization
	r = g_defr;
	minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	oe_del = _mm_set1_epi8(_o_del + _e_del);
	e_del = _mm_set1_epi8(_e_del);
	oe_ins = _mm_set1_epi8(_o_ins + _e_ins);
	e_ins = _mm_set1_epi8(_e_ins);
	shift = _mm_set1_epi8(q->shift);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, cmp, imax;
		__m128i e, h, t, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
		for (j = 0; LIKELY(j < slen); ++j) {
			/* SW cells are computed in the following order:
			 *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			 *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
			 *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
			 */
			// compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
			h = _mm_adds_epu8(h, _mm_load_si128(S + j));
			h = _mm_subs_epu8(h, shift); // h=H'(i-1,j-1)+S(i,j)
			e = _mm_load_si128(E + j); // e=E'(i,j)
			h = _mm_max_epu8(h, e);
			h = _mm_max_epu8(h, f); // h=H'(i,j)
			max = _mm_max_epu8(max, h); // set max
			_mm_store_si128(H1 + j, h); // save to H'(i,j)
			// now compute E'(i+1,j)
			e = _mm_subs_epu8(e, e_del); // e=E'(i,j) - e_del
			t = _mm_subs_epu8(h, oe_del); // h=H'(i,j) - o_del - e_del
			e = _mm_max_epu8(e, t); // e=E'(i+1,j)
			_mm_store_si128(E + j, e); // save to E'(i+1,j)
			// now compute F'(i,j+1)
			f = _mm_subs_epu8(f, e_ins);
			t = _mm_subs_epu8(h, oe_ins); // h=H'(i,j) - o_ins - e_ins
			f = _mm_max_epu8(f, t);
			// get H'(i-1,j) and prepare for the next j
			h = _mm_load_si128(H0 + j); // h=H'(i-1,j)
		}
		// NB: we do not need to set E(i,j) as we disallow adjecent insertion and then deletion
		for (k = 0; LIKELY(k < 16); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
			f = _mm_slli_si128(f, 1);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epu8(h, f); // h=H'(i,j)
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu8(h, oe_ins);
				f = _mm_subs_epu8(f, e_ins);
				cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(f, h), zero));
				if (UNLIKELY(cmp == 0xffff)) goto end_loop16;
			}
		}
end_loop16:
		//int k;for (k=0;k<16;++k)printf("%d ", ((uint8_t*)&max)[k]);printf("\n");
		__max_16(imax, max); // imax is the maximum number in max
		if (imax >= minsc) { // write the b array; this condition adds branching unfornately
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) { // then append
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*)realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i; // te is the end position on the target
			for (j = 0; LIKELY(j < slen); ++j) // keep the H1 vector
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax + q->shift >= 255 || gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S; // swap H0 and H1
	}
	r.score = gmax + q->shift < 255? gmax : 255;
	r.te = te;
	if (r.score != 255) { // get a->qe, the end of query match; find the 2nd best score
		int max = -1, tmp, low, high, qlen = slen * 16;
		uint8_t *t = (uint8_t*)Hmax;
		for (i = 0; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 16 + i % 16 * slen;
			else if ((int)*t == max && (tmp = i / 16 + i % 16 * slen) < r.qe) r.qe = tmp; 
		//printf("%d,%d\n", max, gmax);
		if (b) {
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	free(b);
	return r;
}

kswr_t kswge_i16(kswq_t *q, int tlen, const uint8_t *target, int _o_del, int _e_del, int _o_ins, int _e_ins, int xtra) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i zero, oe_del, e_del, oe_ins, e_ins, *H0, *H1, *E, *Hmax;
	kswr_t r;

#define __max_8(ret, xx) do { \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    	(ret) = _mm_extract_epi16((xx), 0); \
	} while (0)

	// initialization
	r = g_defr;
	minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	oe_del = _mm_set1_epi16(_o_del + _e_del);
	e_del = _mm_set1_epi16(_e_del);
	oe_ins = _mm_set1_epi16(_o_ins + _e_ins);
	e_ins = _mm_set1_epi16(_e_ins);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
		__m128i e, t, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 2);
		for (j = 0; LIKELY(j < slen); ++j) {
			h = _mm_adds_epi16(h, *S++);
			e = _mm_load_si128(E + j);
			h = _mm_max_epi16(h, e);
			h = _mm_max_epi16(h, f);
			max = _mm_max_epi16(max, h);
			_mm_store_si128(H1 + j, h);
			e = _mm_subs_epu16(e, e_del);
			t = _mm_subs_epu16(h, oe_del);
			e = _mm_max_epi16(e, t);
			_mm_store_si128(E + j, e);
			f = _mm_subs_epu16(f, e_ins);
			t = _mm_subs_epu16(h, oe_ins);
			f = _mm_max_epi16(f, t);
			h = _mm_load_si128(H0 + j);
		}
		for (k = 0; LIKELY(k < 16); ++k) {
			f = _mm_slli_si128(f, 2);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epi16(h, f);
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu16(h, oe_ins);
				f = _mm_subs_epu16(f, e_ins);
				if(UNLIKELY(!_mm_movemask_epi8(_mm_cmpgt_epi16(f, h)))) goto end_loop8;
			}
		}
end_loop8:
		__max_8(imax, max);
		if (imax >= minsc) {
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*)realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i;
			for (j = 0; LIKELY(j < slen); ++j)
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S;
	}
	r.score = gmax; r.te = te;
	{
		int max = -1, tmp, low, high, qlen = slen * 8;
		uint16_t *t = (uint16_t*)Hmax;
		for (i = 0, r.qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 8 + i % 8 * slen;
			else if ((int)*t == max && (tmp = i / 8 + i % 8 * slen) < r.qe) r.qe = tmp; 
		if (b) {
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	free(b);
	return r;
}

static inline void revseq(int l, uint8_t *s)
{
	int i, t;
	for (i = 0; i < l>>1; ++i)
		t = s[i], s[i] = s[l - 1 - i], s[l - 1 - i] = t;
}

kswr_t kswge_align2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int xtra, kswq_t **qry)
{
	int size;
	kswq_t *q;
	kswr_t r, rr;
	kswr_t (*func)(kswq_t*, int, const uint8_t*, int, int, int, int, int);

	q = (qry && *qry)? *qry : kswge_qinit((xtra&KSW_XBYTE)? 1 : 2, qlen, query, m, mat);
	if (qry && *qry == 0) *qry = q;
	func = q->size == 2? kswge_i16 : kswge_u8;
	size = q->size;
	r = func(q, tlen, target, o_del, e_del, o_ins, e_ins, xtra);
	if (qry == 0) free(q);
	if ((xtra&KSW_XSTART) == 0 || ((xtra&KSW_XSUBO) && r.score < (xtra&0xffff))) return r;
	revseq(r.qe + 1, query); revseq(r.te + 1, target); // +1 because qe/te points to the exact end, not the position after the end
	q = kswge_qinit(size, r.qe + 1, query, m, mat);
	rr = func(q, tlen, target, o_del, e_del, o_ins, e_ins, KSW_XSTOP | r.score);
	revseq(r.qe + 1, query); revseq(r.te + 1, target);
	free(q);
	if (r.score == rr.score)
		r.tb = r.te - rr.te, r.qb = r.qe - rr.qe;
	return r;
}

kswr_t kswge_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra, kswq_t **qry)
{
	return kswge_align2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, xtra, qry);
}

/********************
 *** SW extension ***
 ********************/

typedef struct {
	int32_t h, e;
} eh_t;

int kswge_extend2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
	eh_t *eh; // score array
	int8_t *qp; // query profile
	int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ins, max_del, max_ie, gscore, max_off;
	assert(h0 > 0);
	// allocate memory
	qp = (int8_t*)malloc(qlen * m);
	eh = (eh_t*)calloc(qlen + 1, 8);
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	// fill the first row
	eh[0].h = h0; eh[1].h = h0 > oe_ins? h0 - oe_ins : 0;
	for (j = 2; j <= qlen && eh[j-1].h > e_ins; ++j)
		eh[j].h = eh[j-1].h - e_ins;
	// adjust $w if it is too large
	k = m * m;
	for (i = 0, max = 0; i < k; ++i) // get the max score
		max = max > mat[i]? max : mat[i];
	max_ins = (int)((double)(qlen * max + end_bonus - o_ins) / e_ins + 1.);
	max_ins = max_ins > 1? max_ins : 1;
	w = w < max_ins? w : max_ins;
	max_del = (int)((double)(qlen * max + end_bonus - o_del) / e_del + 1.);
	max_del = max_del > 1? max_del : 1;
	w = w < max_del? w : max_del; // TODO: is this necessary?
	// DP loop
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
	max_off = 0;
	beg = 0, end = qlen;
	for (i = 0; LIKELY(i < tlen); ++i) {
		int t, f = 0, h1, m = 0, mj = -1;
		int8_t *q = &qp[target[i] * qlen];
		// apply the band and the constraint (if provided)
		if (beg < i - w) beg = i - w;
		if (end > i + w + 1) end = i + w + 1;
		if (end > qlen) end = qlen;
		// compute the first column
		if (beg == 0) {
			h1 = h0 - (o_del + e_del * (i + 1));
			if (h1 < 0) h1 = 0;
		} else h1 = 0;
		for (j = beg; LIKELY(j < end); ++j) {
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Similar to SSE2-SW, cells are computed in the following order:
			//   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
			eh_t *p = &eh[j];
			int h, M = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
			p->h = h1;          // set H(i,j-1) for the next row
			M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
			h = M > e? M : e;   // e and f are guaranteed to be non-negative, so h>=0 even if M<0
			h = h > f? h : f;
			h1 = h;             // save H(i,j) to h1 for the next column
			mj = m > h? mj : j; // record the position where max score is achieved
			m = m > h? m : h;   // m is stored at eh[mj+1]
			t = M - oe_del;
			t = t > 0? t : 0;
			e -= e_del;
			e = e > t? e : t;   // computed E(i+1,j)
			p->e = e;           // save E(i+1,j) for the next row
			t = M - oe_ins;
			t = t > 0? t : 0;
			f -= e_ins;
			f = f > t? f : t;   // computed F(i,j+1)
		}
		eh[end].h = h1; eh[end].e = 0;
		if (j == qlen) {
			max_ie = gscore > h1? max_ie : i;
			gscore = gscore > h1? gscore : h1;
		}
		if (m == 0) break;
		if (m > max) {
			max = m, max_i = i, max_j = mj;
			max_off = max_off > abs(mj - i)? max_off : abs(mj - i);
		} else if (zdrop > 0) {
			if (i - max_i > mj - max_j) {
				if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
			} else {
				if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
			}
		}
		// update beg and end for the next round
		for (j = beg; LIKELY(j < end) && eh[j].h == 0 && eh[j].e == 0; ++j);
		beg = j;
		for (j = end; LIKELY(j >= beg) && eh[j].h == 0 && eh[j].e == 0; --j);
		end = j + 2 < qlen? j + 2 : qlen;
		//beg = 0; end = qlen; // uncomment this line for debugging
	}
	free(eh); free(qp);
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	if (_gtle) *_gtle = max_ie + 1;
	if (_gscore) *_gscore = gscore;
	if (_max_off) *_max_off = max_off;
	return max;
}

int kswge_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off)
{
	return kswge_extend2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w, end_bonus, zdrop, h0, qle, tle, gtle, gscore, max_off);
}

/********************
 * Global alignment *
 ********************/

#define MINUS_INF -0x40000000

static inline uint32_t *push_cigar(int *n_cigar, int *m_cigar, uint32_t *cigar, int op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = (uint32_t*)realloc(cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

int kswge_global2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int *n_cigar_, uint32_t **cigar_)
{
	eh_t *eh;
	int8_t *qp; // query profile
	int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, score, n_col;
	uint8_t *z; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be a little more complex
	if (n_cigar_) *n_cigar_ = 0;
	// allocate memory
	n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
	z = n_cigar_ && cigar_? (uint8_t*)malloc((long)n_col * tlen) : 0;
	qp = (int8_t*)malloc(qlen * m);
	eh = (eh_t*)calloc(qlen + 1, 8);
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	// fill the first row
	eh[0].h = 0; eh[0].e = MINUS_INF;
	for (j = 1; j <= qlen && j <= w; ++j)
		eh[j].h = -(o_ins + e_ins * j), eh[j].e = MINUS_INF;
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = MINUS_INF; // everything is -inf outside the band
	// DP loop
	for (i = 0; LIKELY(i < tlen); ++i) { // target sequence is in the outer loop
		int32_t f = MINUS_INF, h1, beg, end, t;
		int8_t *q = &qp[target[i] * qlen];
		beg = i > w? i - w : 0;
		end = i + w + 1 < qlen? i + w + 1 : qlen; // only loop through [beg,end) of the query sequence
		h1 = beg == 0? -(o_del + e_del * (i + 1)) : MINUS_INF;
		if (n_cigar_ && cigar_) {
			uint8_t *zi = &z[(long)i * n_col];
			for (j = beg; LIKELY(j < end); ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   M(i,j)   = H(i-1,j-1) + S(i,j)
				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
				// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
				// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
				// In practice, this should happen very rarely given a reasonable scoring system.
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				m += q[j];
				d = m >= e? 0 : 1;
				h = m >= e? m : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				t = m - oe_del;
				e -= e_del;
				d |= e > t? 1<<2 : 0;
				e  = e > t? e    : t;
				p->e = e;
				t = m - oe_ins;
				f -= e_ins;
				d |= f > t? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > t? f    : t;
				zi[j - beg] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			for (j = beg; LIKELY(j < end); ++j) {
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				p->h = h1;
				m += q[j];
				h = m >= e? m : e;
				h = h >= f? h : f;
				h1 = h;
				t = m - oe_del;
				e -= e_del;
				e  = e > t? e : t;
				p->e = e;
				t = m - oe_ins;
				f -= e_ins;
				f  = f > t? f : t;
			}
		}
		eh[end].h = h1; eh[end].e = MINUS_INF;
	}
	score = eh[qlen].h;
	if (n_cigar_ && cigar_) { // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0;
		uint32_t *cigar = 0, tmp;
		i = tlen - 1; k = (i + w + 1 < qlen? i + w + 1 : qlen) - 1; // (i,k) points to the last cell
		while (i >= 0 && k >= 0) {
			which = z[(long)i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			if (which == 0)      cigar = push_cigar(&n_cigar, &m_cigar, cigar, 0, 1), --i, --k;
			else if (which == 1) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, 1), --i;
			else                 cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, 1), --k;
		}
		if (i >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, i + 1);
		if (k >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, k + 1);
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;
	}
	free(eh); free(qp); free(z);
	return score;
}

int kswge_global(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *n_cigar_, uint32_t **cigar_)
{
	return kswge_global2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w, n_cigar_, cigar_);
}

int8_t* kswge_gen_smat(const int match_score, const int mismatch_score){
    int8_t* score_mat = (int8_t*)malloc(25*sizeof(int8_t));
    int i = 0, k = 0, j = 0;
    for(; i < 4; ++i){
        for(j = 0; j < 4; ++j){
            score_mat[k++] = i == j ? match_score : -mismatch_score;
        }
        score_mat[k++] = 0; // ambiguous base: no penalty
    }
    for(j = 0; j < 5; ++j) score_mat[k++] = 0;
    return score_mat;
}

uint32_t* kswge_add_cigar(uint32_t* new_cigar, int32_t* p, int32_t* s, uint32_t length, char op){
    if((*p) >= (*s)){
        ++(*s);
        kroundup32(*s);
        new_cigar = (uint32_t*)realloc(new_cigar, (*s)*sizeof(uint32_t));
    }
    new_cigar[(*p) ++] = kswge_gen_cigar(length, op);
    return new_cigar;
}

// choice 0: current not M, 1: current match, 2: current mismatch
uint32_t* store_previous_m(int8_t choice, uint32_t* length_m, uint32_t* length_x, int32_t* p, int32_t* s, uint32_t* new_cigar){
    if((*length_m) && (choice == 2 || !choice)){
        new_cigar = kswge_add_cigar (new_cigar, p, s, (*length_m), '=');
        (*length_m) = 0;
    }else if((*length_x) && (choice == 1 || !choice)){
        new_cigar = kswge_add_cigar (new_cigar, p, s, (*length_x), 'X');
        (*length_x) = 0;
    }
    return new_cigar;
}

void kswge_mark_mismatch(kswr_t* ret, uint8_t* target, int tlen, uint8_t* query, int qlen){
    int32_t mmlen = 0, p = 0, i, oplen, j, s = ret->ncigar + 2, rlens = 0, qlens = 0;
    uint32_t *new_cigar = (uint32_t*)malloc(s*sizeof(uint32_t));
    uint32_t oplen_m = 0,  oplen_x = 0;
    char opchr;
    ret->nvar = ret->nmm = ret->nins = ret->ndel = ret->lsc = ret->tsc = ret->mm = 0;
    target += ret->tb;
    query += ret->qb;
    rlens += ret->tb;
    qlens += ret->qb;
    if(ret->qb > 0){
        new_cigar[p++] = kswge_gen_cigar(ret->qb, 'S');
        ret->lsc = ret->qb;
    }
    bool valid = true;
    for(i = 0; i < ret->ncigar; ++i){
        opchr = kswge_cigar_opchr(ret->cigar[i]);
        oplen = kswge_cigar_oplen(ret->cigar[i]);
        if(opchr == 'M') {
            rlens += oplen;
            qlens += oplen;
            if(rlens > tlen){
                valid = false;
                break;
            }
            for(j = 0; j < oplen; ++j){
                if(*target != *query){
                    ++mmlen;
                    // the previous is match; however the current one is mismatche
                    new_cigar = store_previous_m(2, &oplen_m, &oplen_x, &p, &s, new_cigar);
                    ++oplen_x;
                }else{
                    // the previous is mismatch; however the current one is matche
                    if(oplen_m == 0 && oplen_x != 0){ ++ ret->nmm; ++ ret->nvar; }
                    new_cigar = store_previous_m (1, &oplen_m, &oplen_x, &p, &s, new_cigar);
                    ++oplen_m;
                }
                ++target;
                ++query;
            }
        }else if (opchr == 'I'){
            qlens += oplen;
            query += oplen;
            mmlen += oplen;
            new_cigar = store_previous_m(0, &oplen_m, &oplen_x, &p, &s, new_cigar);
            new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'I');
            ++ret->nins;
            ++ret->nvar;
        }else if(opchr == 'D'){
            rlens += oplen;
            if(rlens > tlen){
                valid = false;
                break;
            }
            target += oplen;
            mmlen += oplen;
            new_cigar = store_previous_m(0, &oplen_m, &oplen_x, &p, &s, new_cigar);
            new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'D');
            ++ret->ndel;
            ++ret->nvar;
        }
    }
    if(!valid){// template length invalid, alignment BAD
        free(new_cigar);
        free(ret->cigar);
        ret->cigar = NULL;
        ret->ncigar = 0;
        ret->score = 0;
        ret->tb = ret->te = ret->qb = ret->qe = -1;
        return;
    }
    if(oplen_m == 0 && oplen_x != 0){ ++ret->nmm; ++ret->nvar; }
    new_cigar = store_previous_m(0, &oplen_m, &oplen_x, &p, &s, new_cigar);

    oplen = qlen - ret->qe - 1;
    if(oplen > 0){
        new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'S');
        ret->tsc = oplen;
        qlens += oplen;
    }
    if(qlens != qlen){// query length invalid, alignment BAD
        free(new_cigar);
        free(ret->cigar);
        ret->cigar = NULL;
        ret->ncigar = 0;
        ret->score = 0;
        ret->tb = ret->te = ret->qb = ret->qe = -1;
    }else{
        ret->ncigar = p;
        free(ret->cigar);
        ret->cigar = new_cigar;
        ret->mm = mmlen;
    }
}

void kswge_ret_output(FILE* f, const kswr_t* a, const char* rseq, const char* qseq, const unsigned char* table, const char* title, int wlen){
    fprintf(f, "%s\n", title);
    fprintf(f, "top1 alignment score: %d\n", a->score);
    fprintf(f, "top2 alignment score: %d\n", a->score2);
    fprintf(f, "ref match range[0 based inclusive]: [%d, %d]\n", a->tb, a->te);
    fprintf(f, "qry match range[0 based inclusive]: [%d, %d]\n", a->qb, a->qe);
    if(a->cigar){
        int32_t c = 0, left = 0, e = 0, qb = a->tb, pb = a->qb;
        uint32_t i;
        while(e < a->ncigar || left > 0){
            int count = 0;
            int32_t q = qb;
            int32_t p = pb;
            // output ref seq
            fprintf(f, "Ref:%8d    ", q + 1);
            for(c = e; c < a->ncigar; ++c){
                char letter = kswge_cigar_opchr(a->cigar[c]);
                uint32_t length = kswge_cigar_oplen(a->cigar[c]);
                uint32_t l = (count == 0 && left > 0) ? left: length;
                for(i = 0; i < l; ++i){
                    if(letter == 'I' || letter == 'S'){
                        fprintf(f, "-");
                        if(letter == 'S' && c == 0){
                            --pb;;
                            --p;
                        };
                    }
                    else{
                        fprintf(f, "%c", *(rseq + q));
                        ++q;
                    }
                    ++count;
                    if(count == wlen) goto step2;
                }
            }
step2:// output match markers
            fprintf(f, "    %d\n                ", q);
            q = qb;
            count = 0;
            for(c = e; c < a->ncigar; ++c){
                char letter = kswge_cigar_opchr(a->cigar[c]);
                uint32_t length = kswge_cigar_oplen(a->cigar[c]);
                uint32_t l = (count == 0 && left > 0) ? left: length;
                for(i = 0; i < l; ++i){
                    if(letter == 'M' || letter == 'X' || letter == '='){
                        if(table[(int)*(rseq + q)] == table[(int)*(qseq + p)])fprintf(f, "|");
                        else fprintf(f, "*");
                        ++q;
                        ++p;
                    }else{
                        fprintf(f, "*");
                        if(letter == 'I' || letter == 'S') ++p;
                        else ++q;
                    }
                    ++count;
                    if(count == wlen){
                        qb = q;
                        goto step3;
                    }
                }
            }
step3:// output query
            p = pb;
            fprintf(f, "\nQry:%8d    ", p + 1);
            count = 0;
            for(c = e; c < a->ncigar; ++c){
                char letter = kswge_cigar_opchr(a->cigar[c]);
                uint32_t length = kswge_cigar_oplen(a->cigar[c]);
                uint32_t l = (count == 0 && left > 0) ? left: length;
                for(i = 0; i < l; ++i){
                    if(letter == 'D') fprintf(f, "-");
                    else{
                        fprintf(f, "%c", *(qseq + p));
                        ++p;
                    }
                    ++count;
                    if(count == wlen){
                        pb = p;
                        left = l - i - 1;
                        e = (left == 0) ? (c + 1) : c;
                        goto end;
                    }
                }
            }
            e = c;
            left = 0;
end:
            fprintf(f, "    %d\n\n", p);
        }
    }
}

kswr_t* kswge_semi_global(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape){
    return kswge_semi_global2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape);
}

kswr_t* kswge_semi_global2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins){
    kswr_t rgl = g_defr;
    rgl = kswge_align2(qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, KSW_XSTART, NULL);
    kswr_t* ret = (kswr_t*)calloc(1, sizeof(kswr_t));
    if(rgl.score <= 0) return ret; // unmapped
    // compute the band-width
    int w = MAX(SW_BW, qlen + tlen);
    // do SW
    int score = kswge_global2(rgl.qe-rgl.qb+1, query+rgl.qb, rgl.te-rgl.tb+1, target+rgl.tb, 5, mat, o_del, e_del, o_ins, e_ins, w, &ret->ncigar, &ret->cigar);
    if(score <= 0){ // unmapped because optimal score not reached
        if(ret->cigar) free(ret->cigar);
        ret->cigar = NULL;
        ret->ncigar = 0;
        ret->qb = ret->qe = ret->tb = ret->te = -1;
    }else{
        ret->score = score;
        ret->qb = rgl.qb;
        ret->tb = rgl.tb;
        ret->qe = rgl.qe;
        ret->te = rgl.te;
    }
    return ret;
}

void kswge_left_align(kswr_t* ret, uint8_t* qints, uint8_t* rints, int beg, int end){
    if(ret->ndel == 0 && ret->nins == 0) return;
    if(beg < 0 || end < 0){ beg = 0; end = ret->te; }
    char opchr; int oplen;
    int rl = ret->te;
    int ql = ret->qe;
    int s = ret->ncigar + 2;
    uint32_t* new_cigar = (uint32_t*)malloc(s*sizeof(uint32_t));
    int p = 0;
    bool keep = true;
    for(int32_t i = ret->ncigar - 1; i >= 0; --i){
        opchr = kswge_cigar_opchr(ret->cigar[i]);
        oplen = kswge_cigar_oplen(ret->cigar[i]);
        if(opchr == '=' || opchr == 'X'){
            ql -= oplen;
            rl -= oplen;
            new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
        }else if(opchr == 'I'){
            if(rl < beg || rl > end){
                new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
                ql -= oplen;
                continue;
            }
            keep = true;
            if(i > 0 && kswge_cigar_opchr(ret->cigar[i-1]) == '='){ // *M*I
                int opnen = kswge_cigar_oplen(ret->cigar[i-1]);
                if(opnen >= oplen){
                    int rpn = 1;
                    for(int c = 0; c < oplen; ++c) if(qints[ql-c] != qints[ql-oplen-c]){ --rpn; break; }
                    if(rpn == 1){
                        ++rpn;
                        bool rpt = true;
                        while(rpt && oplen*rpn <= opnen){
                            for(int c = 0; c < oplen; ++c){
                                if(qints[ql-c] != qints[ql-oplen*rpn-c]){
                                    rpt = false;
                                    break;
                                }
                            }
                            if(rpt) ++rpn;
                        }
                        --rpn;
                        keep = false;
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpn*oplen, '=');
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'I');
                        ret->smask |= KSW_FLALN;
                        rl -= rpn*oplen;
                        ql -= (rpn+1)*oplen;
                        if(opnen - rpn*oplen >= 0) ret->cigar[i-1] = kswge_gen_cigar(opnen - rpn*oplen, '=');
                        else --i;
                    }else{// partial dup
                        int rpl = 0;
                        for(int c = 0; c < oplen; ++c){
                            if(qints[ql-c] == qints[ql-oplen-c]) ++rpl;
                            else break;
                        }
                        if(rpl){
                            keep = false;
                            new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpl, '=');
                            new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'I');
                            ret->smask |= KSW_FLALN;
                            rl -= rpl;
                            ql -= (rpl + oplen);
                        }
                        if(opnen >= rpl) ret->cigar[i-1] = kswge_gen_cigar(opnen-rpl, '=');
                        else --i;
                    }
                }else{
                    int rpl = 0;
                    for(int c = 0; c < opnen; ++c){
                        if(qints[ql-c] == qints[ql-oplen-c]) ++rpl;
                        else break;
                    }
                    if(rpl){
                        keep = false;
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpl, '=');
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'I');
                        ret->smask |= KSW_FLALN;
                        rl -= rpl;
                        ql -= (rpl + oplen);
                    }
                    if(opnen >= rpl) ret->cigar[i-1] = kswge_gen_cigar(opnen-rpl, '=');
                    else --i;
                }
            }
            if(keep){
                new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
                ql -= oplen;
            }
        }else if(opchr == 'D'){
            if(rl - oplen < beg || rl - oplen > end){
                new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
                rl -= oplen;
                continue;
            }
            keep = true;
            if(i > 0 && kswge_cigar_opchr(ret->cigar[i-1]) == '='){ // *M*D
                int opnen = kswge_cigar_oplen(ret->cigar[i-1]);
                if(opnen >= oplen){
                    int rpn = 1;
                    for(int c = 0; c < oplen; ++c) if(rints[rl-c] != rints[rl-oplen-c]){ --rpn; break; }
                    if(rpn == 1){
                        ++rpn;
                        bool rpt = true;
                        while(rpt && oplen*rpn < opnen){
                            for(int c = 0; c < oplen; ++c){
                                if(rints[rl-c] != rints[rl-oplen*rpn-c]){
                                    rpt = false;
                                    break;
                                }
                            }
                            if(rpt) ++rpn;
                        }
                        --rpn;
                        keep = false;
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpn*oplen, '=');
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'D');
                        ret->smask |= KSW_FLALN;
                        ql -= rpn * oplen;
                        rl -= (rpn+1)*oplen;
                        if(opnen - rpn*oplen >= 0) ret->cigar[i-1] = kswge_gen_cigar(opnen - rpn*oplen, '=');
                        else --i;
                    }else{ // partial dup
                        int rpl = 0;
                        for(int c = 0; c < oplen; ++c){
                            if(rints[rl-c] == rints[rl-oplen-c]) ++rpl;
                            else break;
                        }
                        if(rpl){ // pritial dup found
                            keep = false;
                            new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpl, '=');
                            new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'D');
                            ret->smask |= KSW_FLALN;
                            ql -= rpl;
                            rl -= (rpl + oplen);;
                            if(opnen >= rpl) ret->cigar[i-1] = kswge_gen_cigar(opnen - rpl, '=');
                            else --i;
                        }
                    }
                }else{
                    int rpl = 0;
                    for(int c = 0; c < opnen; ++c){
                        if(rints[rl-c] == rints[rl-oplen-c]) ++rpl;
                        else break;
                    }
                    if(rpl){ // pritial dup found
                        keep = false;
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpl, '=');
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'D');
                        ret->smask |= KSW_FLALN;
                        ql -= rpl;
                        rl -= (rpl+oplen);
                        if(opnen >= rpl) ret->cigar[i-1] = kswge_gen_cigar(opnen - rpl, '=');
                        else --i;
                    }
                }
            }
            if(keep){
                new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
                rl -= oplen;
            }
        }else if(opchr == 'S'){
            new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
        }
    }
    // reverse cigar to normal order
    int bpos = 0, epos = p-1;
    while(bpos < epos){
        new_cigar[bpos] = new_cigar[bpos] ^ new_cigar[epos];
        new_cigar[epos] = new_cigar[bpos] ^ new_cigar[epos];
        new_cigar[bpos] = new_cigar[bpos] ^ new_cigar[epos];
        ++bpos;
        --epos;
    }
    // merge adjacent operations if possible
    for(int k = 1; k < p; ++k){
        if(kswge_cigar_opchr(new_cigar[k]) == kswge_cigar_opchr(new_cigar[k-1])){
            new_cigar[k] += new_cigar[k-1] >> BAM_CIGAR_SHIFT << BAM_CIGAR_SHIFT, new_cigar[k-1] &= 0xf;
        }
    }
    // kill zero length operations
    int k = 0, i = 0;
    for(k = 0, i = 0; k < p; ++k){
        if(new_cigar[k] >> BAM_CIGAR_SHIFT) new_cigar[i++] = new_cigar[k];
    }
    // update var count
    ret->nvar = ret->nmm = ret->nins = ret->ndel = 0;
    for(k = 0; k < i; ++k){
        switch(kswge_cigar_opchr(new_cigar[k])){
            case 'X':
                ++ret->nmm;
                ++ret->nvar;
                break;
            case 'D':
                ++ret->ndel;
                ++ret->nvar;
                break;
            case 'I':
                ++ret->nins;
                ++ret->nvar;
                break;
            default:
                break;
        }
    }
    ret->ncigar = i;
    free(ret->cigar);
    ret->cigar = new_cigar;
}

void kswge_right_align(kswr_t* ret, uint8_t* qints, uint8_t* rints, int beg, int end){
    if(ret->ndel == 0 && ret->nins == 0) return;
    if(beg < 0 || end < 0){ beg = 0; end = ret->te; }
    char opchr; int oplen;
    int rl = ret->tb;
    int ql = ret->qb;
    int s = ret->ncigar + 2;
    uint32_t* new_cigar = (uint32_t*)malloc(s*sizeof(uint32_t));
    int p = 0;
    bool keep = true;
    for(int32_t i = 0; i < ret->ncigar; ++i){
        opchr = kswge_cigar_opchr(ret->cigar[i]);
        oplen = kswge_cigar_oplen(ret->cigar[i]);
        if(opchr == '=' || opchr == 'X'){
            ql += oplen;
            rl += oplen;
            new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
        }else if(opchr == 'I'){
            if(rl < beg || rl > end){
                new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
                ql += oplen;
                continue;
            }
            keep = true;
            if(i < ret->ncigar - 1 && kswge_cigar_opchr(ret->cigar[i+1]) == '='){ // *I*M
                int opnen = kswge_cigar_oplen(ret->cigar[i+1]);
                if(opnen >= oplen){
                    int rpn = 1;
                    for(int c = 0; c < oplen; ++c) if(qints[ql+c] != qints[ql+oplen+c]){ --rpn; break; }
                    if(rpn == 1){
                        ++rpn;
                        bool rpt = true;
                        while(rpt && oplen*rpn <= opnen){
                            for(int c = 0; c < oplen; ++c){
                                if(qints[ql+c] != qints[ql+oplen*rpn+c]){
                                    rpt = false;
                                    break;
                                }
                            }
                            if(rpt) ++rpn;
                        }
                        --rpn;
                        keep = false;
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpn*oplen, '=');
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'I');
                        ret->smask |= KSW_FRALN;
                        rl += rpn * oplen;
                        ql += (rpn+1)*oplen;
                        if(opnen - rpn*oplen >= 0) ret->cigar[i+1] = kswge_gen_cigar(opnen - rpn*oplen, '=');
                        else ++i;
                    }else{// partial dup
                        int rpl = 0;
                        for(int c = 0; c < oplen; ++c){
                            if(qints[ql+c] == qints[ql+oplen+c]) ++rpl;
                            else break;
                        }
                        if(rpl){// partial dup found
                            keep = false;
                            new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpl, '=');
                            new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'I');
                            ret->smask |= KSW_FRALN;
                            rl += rpl;
                            ql += rpl + oplen;
                            if(opnen >= rpl) ret->cigar[i+1] = kswge_gen_cigar(opnen - rpl, '=');
                            else ++i;
                        }
                    }
                }else{
                    int rpl = 0;
                    for(int c = 0; c < opnen; ++c){
                        if(qints[ql+c] == qints[ql+oplen+c]) ++rpl;
                        else break;
                    }
                    if(rpl){// partial dup found
                        keep = false;
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpl, '=');
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'I');
                        ret->smask |= KSW_FRALN;
                        rl += rpl;
                        ql += rpl + oplen;
                    }
                    if(opnen >= rpl) ret->cigar[i+1] = kswge_gen_cigar(opnen-rpl, '=');
                    else ++i;
                }
            }
            if(keep){
                new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
                ql += oplen;
            }
        }else if(opchr == 'D'){
            if(rl < beg || rl > end){
                rl += oplen;
                new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
                continue;
            }
            keep = true;
            if(i < ret->ncigar - 1 && kswge_cigar_opchr(ret->cigar[i+1]) == '='){ // *D*M
                int opnen = kswge_cigar_oplen(ret->cigar[i+1]);
                if(opnen >= oplen){
                    int rpn = 1;
                    for(int c = 0; c < oplen; ++c) if(rints[rl+c] != rints[rl+oplen+c]){ --rpn; break; }
                    if(rpn == 1){
                        ++rpn;
                        bool rpt = true;
                        while(rpt && oplen*rpn < opnen){
                            for(int c = 0; c < oplen; ++c){
                                if(rints[rl+c] != rints[rl+oplen*rpn+c]){
                                    rpt = false;
                                    break;
                                }
                            }
                            if(rpt) ++rpn;
                        }
                        --rpn;
                        keep = false;
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpn*oplen, '=');
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'D');
                        ret->smask |= KSW_FRALN;
                        ql += rpn*oplen;
                        rl += (rpn+1)*oplen;
                        if(opnen - rpn*oplen >= 0) ret->cigar[i+1] = kswge_gen_cigar(opnen - rpn*oplen, '=');
                        else ++i;
                    }else{ // partial dup
                        int rpl = 0;
                        for(int c = 0; c < oplen; ++c){
                            if(rints[rl+c] == rints[rl+oplen+c]) ++rpl;
                            else break;
                        }
                        if(rpl){ // pritial dup found
                            keep = false;
                            new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpl, '=');
                            new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'D');
                            ret->smask |= KSW_FRALN;
                            ql += rpl;
                            rl += rpl+oplen;
                            if(opnen >= rpl) ret->cigar[i+1] = kswge_gen_cigar(opnen - rpl, '=');
                            else ++i;
                        }
                    }
                }else{
                    int rpl = 0;
                    for(int c = 0; c < opnen; ++c){
                        if(rints[rl+c] == rints[rl+oplen+c]) ++rpl;
                        else break;
                    }
                    if(rpl){ // pritial dup found
                        keep = false;
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, rpl, '=');
                        new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, 'D');
                        ret->smask |= KSW_FRALN;
                        ql += rpl;
                        rl += rpl + oplen;
                        if(opnen >= rpl) ret->cigar[i+1] = kswge_gen_cigar(opnen - rpl, '=');
                        else ++i;
                    }
                }
            }
            if(keep){
                new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
                rl += oplen;
            }
        }else if(opchr == 'S'){
            new_cigar = kswge_add_cigar(new_cigar, &p, &s, oplen, opchr);
        }
    }
    // merge adjacent operations if possible
    for(int k = 1; k < p; ++k){
        if(kswge_cigar_opchr(new_cigar[k]) == kswge_cigar_opchr(new_cigar[k-1])){
            new_cigar[k] += new_cigar[k-1] >> BAM_CIGAR_SHIFT << BAM_CIGAR_SHIFT, new_cigar[k-1] &= 0xf;
        }
    }
    // kill zero length operations
    int k = 0, i = 0;
    for(k = 0, i = 0; k < p; ++k){
        if(new_cigar[k] >> BAM_CIGAR_SHIFT) new_cigar[i++] = new_cigar[k];
    }
    // update var count
    ret->nvar = ret->nmm = ret->nins = ret->ndel = 0;
    for(k = 0; k < i; ++k){
        switch(kswge_cigar_opchr(new_cigar[k])){
            case 'X':
                ++ret->nmm;
                ++ret->nvar;
                break;
            case 'D':
                ++ret->ndel;
                ++ret->nvar;
                break;
            case 'I':
                ++ret->nins;
                ++ret->nvar;
                break;
            default:
                break;
        }
    }
    ret->ncigar = i;
    free(ret->cigar);
    ret->cigar = new_cigar;
}

int kswge_indel_compatible(kswr_t* aln1, kswr_t* aln2, const char* qseq1, const char* qseq2){
    // r1...r2 || r2...r1 big ins
    if(aln1->cigar == NULL || aln2->cigar == NULL) return -1; // single unmap
    if(aln1->te < aln2->tb || aln1->tb > aln2->te) return 0; // non overlap
    int tm = MAX(aln1->te + 2, aln2->te + 2);
    char** insseq1 = (char**)calloc(tm, sizeof(char*));
    int* inslen1 = (int*)calloc(tm, sizeof(int));
    int* dellen1 = (int*)calloc(tm, sizeof(int));
    int* dislen1 = (int*)calloc(tm, sizeof(int));
    char** disseq1 = (char**)calloc(tm, sizeof(char*));
    char** insseq2 = (char**)calloc(tm, sizeof(char*));
    int* inslen2 = (int*)calloc(tm, sizeof(int));
    int* dellen2 = (int*)calloc(tm, sizeof(int));
    int* dislen2 = (int*)calloc(tm, sizeof(int));
    char** disseq2 = (char**)calloc(tm, sizeof(char*));
    int oplen;
    char opchr;
    int rpos = aln1->tb;
    int qpos = aln1->qb;
    for(int32_t i = 0; i < aln1->ncigar; ++i){
        opchr = kswge_cigar_opchr(aln1->cigar[i]);
        oplen = kswge_cigar_oplen(aln1->cigar[i]);
        switch(opchr){
            case 'D':
                dellen1[rpos] = oplen;
                rpos += oplen;
                break;
            case 'I':
                inslen1[rpos] = oplen;
                insseq1[rpos] = strndup(qseq1 + qpos, oplen);
                qpos += oplen;
                break;
            case 'X':
                if(oplen > 1){// delins
                    dislen1[rpos] = oplen;
                    disseq1[rpos] = strndup(qseq1 + qpos, oplen);
                }
                rpos += oplen;
                qpos += oplen;
                break;
            case '=':
                rpos += oplen;
                qpos += oplen;
                break;
            default:
                break;
        }
    }
    // check q2
    rpos = aln2->tb;
    qpos = aln2->qb;
    for(int32_t i = 0; i < aln2->ncigar; ++i){
        opchr = kswge_cigar_opchr(aln2->cigar[i]);
        oplen = kswge_cigar_oplen(aln2->cigar[i]);
        switch(opchr){
            case 'D':
                dellen2[rpos] = oplen;
                rpos += oplen;
                break;
            case 'I':
                inslen2[rpos] = oplen;
                insseq2[rpos] = strndup(qseq2 + qpos, oplen);
                qpos += oplen;
                break;
            case 'X':
                if(oplen > 1){// delins
                    dislen2[rpos] = oplen;
                    disseq2[rpos] = strndup(qseq2 + qpos, oplen);
                }
                rpos += oplen;
                qpos += oplen;
                break;
            case '=':
                rpos += oplen;
                qpos += oplen;
                break;
            default:
                break;
        }
    }
    // consistent check
    bool compatible = true;
    for(int i = 0; i < tm; ++i){
        if(dellen1[i] != dellen2[i]){
            compatible = false;
            break;
        }
        if(inslen1[i] != inslen2[i]){
            compatible = false;
            break;
        }else if(inslen1[i] == inslen2[i] && inslen1[i]){
            if(strcmp(insseq1[i], insseq2[i])){
                compatible = false;
                break;
            }
        }
        if(dislen1[i] != dislen2[i]){
            compatible = false;
            break;
        }else if(dislen1[i] == dislen2[i] && dislen1[i]){
            if(strcmp(disseq1[i], disseq2[i])){
                compatible = false;
                break;
            }
        }
    }
    for(int i = 0; i < tm; ++i){
        if(insseq1[i]) free(insseq1[i]);
        if(insseq2[i]) free(insseq2[i]);
        if(disseq1[i]) free(disseq1[i]);
        if(disseq2[i]) free(disseq2[i]);
    }
    free(insseq1);
    free(insseq2);
    free(inslen1);
    free(inslen2);
    free(dellen1);
    free(dellen2);
    free(dislen1);
    free(dislen2);
    free(disseq1);
    free(disseq2);
    if(compatible) return 1;
    else return -1;
}

void kswge_usage(kswge_opt_t* opt, char* arg0){
     fprintf(stderr, "\n");
     fprintf(stderr, "Usage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
     fprintf(stderr, "Options: -m INT match score [%d]\n", opt->match);
     fprintf(stderr, "         -d INT mismatch score [%d]\n", opt->mismatch);
     fprintf(stderr, "         -o INT gap open penalty [%d]\n", opt->gapo);
     fprintf(stderr, "         -e INT gap extend penalty [%d]\n", opt->gape);
     fprintf(stderr, "         -r STR reference sequence\n");
     fprintf(stderr, "         -q STR query sequence\n");
     fprintf(stderr, "         -Q STR query sequence\n");
     fprintf(stderr, "         -1 INT beg left aln pos\n");
     fprintf(stderr, "         -2 INT end left aln pos\n");
     fprintf(stderr, "         -3 INT beg right aln pos\n");
     fprintf(stderr, "         -4 INT end right aln pos\n");
     fprintf(stderr, "         -w INT output line width [%d]\n", opt->wlen);
     fprintf(stderr, "\n");
}

int kswge_main(int argc, char** argv)
{
    kswge_opt_t opt;
    if(argc == 1){
        kswge_usage(&opt, argv[0]);
        return 0;
    }
    int c = 0;
    while((c = getopt(argc, argv, "m:d:o:e:r:q:Q:1:2:3:4:w:h")) >= 0){
        switch(c){
            case 'm': opt.match = atoi(optarg); break;
            case 'd': opt.mismatch = atoi(optarg); break;
            case 'o': opt.gapo = atoi(optarg); break;
            case 'e': opt.gape = atoi(optarg); break;
            case 'r': opt.ref = strdup(optarg); break;
            case 'q': opt.qry = strdup(optarg); break;
            case 'Q': opt.qry2 = strdup(optarg); break;
            case '1': opt.bpos = atoi(optarg); break;
            case '2': opt.epos = atoi(optarg); break;
            case '3': opt.posb = atoi(optarg); break;
            case '4': opt.pose = atoi(optarg); break;
            case 'w': opt.wlen = atoi(optarg); break;
            case 'h': kswge_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    // check
    if(opt.qry == NULL || opt.ref == NULL){
        fprintf(stderr, "-q(qry) and -r(ref) must be provided\n");
        return 1;
    }
    // seq to ints
    int rlen = strlen(opt.ref), qlen = strlen(opt.qry), q2len = 0;
    if(opt.qry2) q2len = strlen(opt.qry2);
    uint8_t* rints = kswge_seq2ints(opt.ref, rlen);
    uint8_t* qints = kswge_seq2ints(opt.qry, qlen);
    uint8_t* q2ints = NULL;
    if(opt.qry2) q2ints = kswge_seq2ints(opt.qry2, q2len);
    // aligment aux
    int xtra = KSW_XSTART;
    int8_t* mat = kswge_gen_smat(opt.match, opt.mismatch);
    // q1
    /** local */
    kswr_t rl = kswge_align(qlen, qints, rlen, rints, 5, mat, opt.gapo, opt.gape, xtra, NULL);
    fprintf(stderr, "local q range of r[%d, %d]\n", rl.qb, rl.qe);
    fprintf(stderr, "local q range of t[%d, %d]\n", rl.tb, rl.te);
    /** global */
    int w = MAX(rlen + qlen, SW_BW);
    kswr_t rg = g_defr;
    int score = kswge_global(qlen, qints, rlen, rints, 5, mat, opt.gapo, opt.gape, w, &rg.ncigar, &rg.cigar);
    /** important */
    rg.tb = rg.qb = 0;
    rg.te = rlen - 1;
    rg.qe = qlen - 1;
    rg.score = score;
    /** important */
    kswge_cigar_output(stderr, rg.cigar, rg.ncigar, "global q before kswge_mark_mismatch :");
    kswge_ret_output(stderr, &rg, opt.ref, opt.qry, kswge_nt2int, "global alignment result of q", opt.wlen);
    kswge_mark_mismatch(&rg, rints, rlen, qints, qlen);
    kswge_cigar_output(stderr, rg.cigar, rg.ncigar, "global q finish kswge_mark_mismatch :");
    kswge_ret_output(stderr, &rg, opt.ref, opt.qry, kswge_nt2int, "global alignment result of q", opt.wlen);
    /** semi global */
    kswr_t* rs = kswge_semi_global(qlen, qints, rlen, rints, 5, mat, opt.gapo, opt.gape);
    kswge_cigar_output(stderr, rs->cigar, rs->ncigar, "semiglobal q before kswge_mark_mismatch :");
    kswge_ret_output(stderr, rs, opt.ref, opt.qry, kswge_nt2int, "semiglobal alignment result of q", opt.wlen);
    kswge_mark_mismatch(rs, rints, rlen, qints, qlen);
    kswge_cigar_output(stderr, rs->cigar, rs->ncigar, "semiglobal q finish kswge_mark_mismatch :");
    kswge_ret_output(stderr, rs, opt.ref, opt.qry, kswge_nt2int, "semiglobal alignment result of q", opt.wlen);
    /** right align */
    kswge_right_align(rs, qints, rints, opt.posb, opt.pose);
    kswge_cigar_output(stderr, rs->cigar, rs->ncigar, "semiglobal q after right align :");
    kswge_ret_output(stderr, rs, opt.ref, opt.qry, kswge_nt2int, "semiglobal alignment result of q", opt.wlen);
    /** left align */
    kswge_left_align(rs, qints, rints, opt.bpos, opt.epos);
    kswge_cigar_output(stderr, rs->cigar, rs->ncigar, "semiglobal q after left align :");
    kswge_ret_output(stderr, rs, opt.ref, opt.qry, kswge_nt2int, "semiglobal alignment result of q", opt.wlen);
    // q2
    if(opt.qry2){
        kswr_t* rs2 = kswge_semi_global(q2len, q2ints, rlen, rints, 5, mat, opt.gapo, opt.gape);
        kswge_cigar_output(stderr, rs2->cigar, rs2->ncigar, "semiglobal Q before kswge_mark_mismatch :");
        kswge_mark_mismatch(rs2, rints, rlen, q2ints, q2len);
        kswge_cigar_output(stderr, rs2->cigar, rs2->ncigar, "semiglobal Q finish kswge_mark_mismatch :");
        /** right align */
        kswge_right_align(rs2, q2ints, rints);
        int comp = kswge_indel_compatible(rs, rs2, opt.qry, opt.qry2);
        if(comp == 0){
            fprintf(stderr, "q and Q not overlap on ref\n");
        }else if(comp == 1){
            fprintf(stderr, "q and Q overlap and has consistent indel\n");
        }else if(comp == -1){
            fprintf(stderr, "q and Q overlap and has inconsistent indel\n");
        }
        if(rs2->cigar) free(rs2->cigar);
        free(rs2);
    }
    // release res
    if(rg.cigar) free(rg.cigar);
    if(rs->cigar) free(rs->cigar);
    free(rs);
    free(opt.qry); free(qints);
    if(opt.qry2){ free(opt.qry2); free(q2ints); }
    free(mat);
    return 0;
}


#ifdef _KSW_MAIN
int main(int argc, char** argv){
    return kswge_main(argc, argv);
}
#endif
