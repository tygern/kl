// Command e6-mu-table computes the complete Kazhdan-Lusztig table of a finite
// simply laced Weyl group of rank at most 8 (types A, D, E) with every leading
// coefficient mu(x,w), in exact integer arithmetic, and writes the E6
// certificate of the manuscript's Section on complements: the maximum of mu
// over all pairs, the pairs attaining it, the histogram of mu >= 2, and the
// fact that every pair with a fully commutative lower endpoint has mu in {0,1}.
//
// Port of research/review_checks/e6_mu_table.cpp (same method, same
// assertions, same output schema).
//
// Method.  Elements are integer matrices in the simple-root basis (column j is
// w(alpha_j)).  The group is enumerated by breadth-first right multiplication,
// Bruhat order is stored as bitsets, and the Kazhdan-Lusztig left recursion
// (Kazhdan-Lusztig 1979, (2.2.c)) is evaluated on extremal pairs only
// (x extremal for w iff x <= w, L(w) subset L(x), R(w) subset R(x)).  Every
// polynomial is recomputed with the right recursion (w = v t, t the largest
// right descent) and the two values must match. Constant terms,
// nonnegativity and the degree bound are asserted for every pair.  mu(x,w) for
// a Bruhat cover x of w is 1 by definition; covers are found as w t over all
// reflections t.
//
// Usage: e6-mu-table -type E -rank 6 -out results/e6-mu-table.json
//
// Hidden debugging flag (differential testing only, no certificate output):
// e6-mu-table -type T -rank n -pairs FILE reads lines "xword wword" (words
// as strings of generator labels, "e" for the identity) and prints one line
// "xword wword P_ascending" per pair with P_{x,w} from the completed table
// ("[]" when x is not below w), then exits without writing any certificate.
// With -pairs absent the program's behaviour and output are unchanged.
//
// Diagram labelling for type E: chain 0-1-...-(n-2), node n-1 attached to
// node 2 (the manuscript's convention); type D: chain 0-...-(n-2), node n-1
// attached to node n-3; type A: chain 0-...-(n-1).
//
// Imports: standard library only.  This command is self-contained and shares
// no code with the other Go engines of the supplement.
package main

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"flag"
	"fmt"
	"math/bits"
	"os"
	"sort"
	"strconv"
	"strings"
)

type mat [8][8]int8 // c[j][k] = coefficient of alpha_k in w(alpha_j)

var (
	n int
	A [8][8]int // Cartan matrix, A[i][j] = <alpha_i^vee, alpha_j>
	M [8][8]int // Coxeter matrix
)

func fail(code int, msg string) {
	fmt.Fprintln(os.Stderr, msg)
	os.Exit(code)
}

func chk(v int) int8 {
	if v > 127 || v < -127 {
		fail(3, "int8 overflow")
	}
	return int8(v)
}

func rightMul(m *mat, i int) mat {
	r := *m
	for j := 0; j < n; j++ {
		if j == i {
			continue
		}
		a := A[i][j]
		if a != 0 {
			for k := 0; k < n; k++ {
				r[j][k] = chk(int(m[j][k]) - a*int(m[i][k]))
			}
		}
	}
	for k := 0; k < n; k++ {
		r[i][k] = chk(-int(m[i][k]))
	}
	return r
}

func leftMul(m *mat, i int) mat {
	r := *m
	for j := 0; j < n; j++ {
		pair := 0
		for k := 0; k < n; k++ {
			pair += A[i][k] * int(m[j][k])
		}
		r[j][i] = chk(int(m[j][i]) - pair)
	}
	return r
}

func colNegative(m *mat, j int) bool {
	for k := 0; k < n; k++ {
		if m[j][k] != 0 {
			return m[j][k] < 0
		}
	}
	return false
}

var (
	N     uint32
	mats  []mat
	slen  []uint8
	rmul  [][8]uint32
	lmul  [][8]uint32
	Lmask []uint8
	Rmask []uint8

	bruhat [][]uint64
)

func leq(x, w uint32) bool { return (bruhat[w][x>>6]>>(x&63))&1 == 1 }

// leqIter is an independent Bruhat test by the lifting property, used to
// cross-check the bitsets.
func leqIter(x, w uint32) bool {
	for {
		if x == w {
			return true
		}
		if slen[x] >= slen[w] {
			return false
		}
		s := bits.TrailingZeros8(Lmask[w])
		w = lmul[w][s]
		if (Lmask[x]>>uint(s))&1 == 1 {
			x = lmul[x][s]
		}
	}
}

// Polynomial table with hash-consing.
var (
	polys     [][]int64
	polyIndex = map[string]uint32{}
	keyBuf    []byte
)

func polyKey(p []int64) []byte {
	keyBuf = keyBuf[:0]
	for _, c := range p {
		keyBuf = binary.LittleEndian.AppendUint64(keyBuf, uint64(c))
	}
	return keyBuf
}

func internPoly(p []int64) uint32 {
	for len(p) > 0 && p[len(p)-1] == 0 {
		p = p[:len(p)-1]
	}
	k := polyKey(p)
	if id, ok := polyIndex[string(k)]; ok {
		return id
	}
	id := uint32(len(polys))
	cp := make([]int64, len(p))
	copy(cp, p)
	polys = append(polys, cp)
	polyIndex[string(k)] = id
	return id
}

type zmu struct {
	z  uint32
	mu int64
}

var (
	ex     [][]uint32 // sorted extremal x for each w
	pid    [][]uint32 // polynomial ids parallel to ex
	mulist [][]zmu    // (z, mu(z,w)) with z < w and mu != 0
)

func climb(x, v uint32) uint32 {
	for {
		dl := Lmask[v] &^ Lmask[x]
		dr := Rmask[v] &^ Rmask[x]
		if dl == 0 && dr == 0 {
			return x
		}
		if dl != 0 {
			x = lmul[x][bits.TrailingZeros8(dl)]
		} else {
			x = rmul[x][bits.TrailingZeros8(dr)]
		}
	}
}

func lookup(x, v uint32) ([]int64, bool) {
	if !leq(x, v) {
		return nil, false
	}
	x = climb(x, v)
	e := ex[v]
	i := sort.Search(len(e), func(i int) bool { return e[i] >= x })
	if i == len(e) || e[i] != x {
		fail(3, "lookup failure")
	}
	return polys[pid[v][i]], true
}

// reducedWord is the canonical reduced word: repeatedly strip the smallest
// right descent.
func reducedWord(x uint32) []int {
	var w []int
	for x != 0 {
		s := bits.TrailingZeros8(Rmask[x])
		w = append(w, s)
		x = rmul[x][s]
	}
	for i, j := 0, len(w)-1; i < j; i, j = i+1, j-1 {
		w[i], w[j] = w[j], w[i]
	}
	return w
}

func wordStr(w []int) string {
	var sb strings.Builder
	for _, a := range w {
		sb.WriteString(strconv.Itoa(a))
	}
	if sb.Len() == 0 {
		return "e"
	}
	return sb.String()
}

func jsonInts[T int | int64](v []T) string {
	var sb strings.Builder
	sb.WriteByte('[')
	for i, a := range v {
		if i > 0 {
			sb.WriteByte(',')
		}
		sb.WriteString(strconv.FormatInt(int64(a), 10))
	}
	sb.WriteByte(']')
	return sb.String()
}

func maskList(m uint8) []int {
	v := []int{}
	for i := 0; i < n; i++ {
		if (m>>uint(i))&1 == 1 {
			v = append(v, i)
		}
	}
	return v
}

func fromWord(w string) uint32 {
	var e uint32
	for _, c := range []byte(w) {
		e = rmul[e][c-'0']
	}
	return e
}

// parseWord reads a word for the -pairs debugging flag: a string of
// generator labels 0..n-1, or "e" for the identity.
func parseWord(w string) uint32 {
	if w == "e" {
		return 0
	}
	for _, c := range []byte(w) {
		if c < '0' || int(c-'0') >= n {
			fail(2, "bad generator label in word "+w)
		}
	}
	return fromWord(w)
}

// dumpPairs implements the hidden -pairs flag: for every line "xword wword"
// of the file it prints "xword wword P" with P = P_{x,w} in ascending powers
// ("[]" if x is not below w).  It is used only for differential testing
// against other implementations and never writes a certificate.
func dumpPairs(path string) {
	f, err := os.Open(path)
	if err != nil {
		fail(2, "cannot open pairs file "+path)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		fields := strings.Fields(sc.Text())
		if len(fields) == 0 || strings.HasPrefix(fields[0], "#") {
			continue
		}
		if len(fields) < 2 {
			fail(2, "pairs line needs two words")
		}
		x, w := parseWord(fields[0]), parseWord(fields[1])
		p, ok := lookup(x, w)
		if !ok {
			p = []int64{}
		}
		fmt.Printf("%s %s %s\n", fields[0], fields[1], jsonInts(p))
	}
	if err := sc.Err(); err != nil {
		fail(2, "reading pairs file: "+err.Error())
	}
}

var failures = 0

func expect(cond bool, what string) {
	if !cond {
		fmt.Fprintf(os.Stderr, "EXPECTATION FAILED: %s\n", what)
		failures++
	}
}

func eqPoly(a, b []int64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func boolStr(b bool) string {
	if b {
		return "true"
	}
	return "false"
}

// guardBound is the overflow guard for stored coefficients and mu values:
// the largest power of two B with N*B*B + 2*B < 2^63, that is
// B = 2^floor((62 - bitlen(N))/2), set once the group order N is known.
// Every polynomial is checked against it before it is interned and every mu
// is read from a checked polynomial, and a recurrence step adds two stored
// polynomials and at most N products mu*coefficient of stored values, so no
// intermediate sum can overflow int64 before guard runs on the result. For
// E6 (N = 51,840) B = 2^23 = 8,388,608.
var guardBound int64

// guard: a polynomial coefficient must stay within the overflow guard.
func guard(p []int64) {
	for _, c := range p {
		if c > guardBound || c < -guardBound {
			fail(3, "coefficient overflow guard")
		}
	}
}

type bigEntry struct {
	x, w uint32
	mu   int64
}

func main() {
	typeFlag := flag.String("type", "", "A, D or E")
	rankFlag := flag.Int("rank", 0, "rank (2..8)")
	outFlag := flag.String("out", "", "output JSON path")
	pairsFlag := flag.String("pairs", "", "debugging: file of \"xword wword\" lines; print P_{x,w} per line and write no certificate")
	flag.Parse()
	if *typeFlag == "" || (*outFlag == "" && *pairsFlag == "") || flag.NArg() != 0 {
		fail(2, "usage: e6-mu-table -type <A|D|E> -rank <n> -out <output.json>")
	}
	typ := (*typeFlag)[0]
	n = *rankFlag
	outfile := *outFlag
	if n < 2 || n > 8 {
		fail(2, "rank must be between 2 and 8")
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				M[i][j] = 1
			} else {
				M[i][j] = 2
			}
		}
	}
	edge := func(i, j int) { M[i][j], M[j][i] = 3, 3 }
	switch typ {
	case 'A':
		for i := 0; i+1 < n; i++ {
			edge(i, i+1)
		}
	case 'D':
		for i := 0; i+1 < n-1; i++ {
			edge(i, i+1)
		}
		edge(n-3, n-1)
	case 'E':
		for i := 0; i+1 < n-1; i++ {
			edge(i, i+1)
		}
		edge(2, n-1)
	default:
		fail(2, "type must be A, D or E")
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				A[i][j] = 2
			} else if M[i][j] == 3 {
				A[i][j] = -1
			}
		}
	}

	// Enumerate the group.
	index := map[mat]uint32{}
	var id mat
	for i := 0; i < n; i++ {
		id[i][i] = 1
	}
	mats = append(mats, id)
	slen = append(slen, 0)
	index[id] = 0
	rmul = append(rmul, [8]uint32{})
	lmul = append(lmul, [8]uint32{})
	for w := 0; w < len(mats); w++ {
		m := mats[w]
		for i := 0; i < n; i++ {
			r := rightMul(&m, i)
			if got, ok := index[r]; !ok {
				id2 := uint32(len(mats))
				mats = append(mats, r)
				slen = append(slen, slen[w]+1)
				index[r] = id2
				rmul = append(rmul, [8]uint32{})
				lmul = append(lmul, [8]uint32{})
				rmul[w][i] = id2
			} else {
				rmul[w][i] = got
			}
		}
	}
	N = uint32(len(mats))
	guardBound = int64(1) << uint((62-bits.Len(uint(N)))/2)
	for w := uint32(0); w < N; w++ {
		for i := 0; i < n; i++ {
			l := leftMul(&mats[w], i)
			got, ok := index[l]
			if !ok {
				fail(3, "left multiplication not closed")
			}
			lmul[w][i] = got
		}
	}
	index = nil
	Lmask = make([]uint8, N)
	Rmask = make([]uint8, N)
	for w := uint32(0); w < N; w++ {
		for i := 0; i < n; i++ {
			if colNegative(&mats[w], i) {
				Rmask[w] |= 1 << uint(i)
			}
			if slen[lmul[w][i]] < slen[w] {
				Lmask[w] |= 1 << uint(i)
			}
		}
	}
	maxLen := int(slen[N-1])
	fmt.Fprintf(os.Stderr, "type %c%d: |W| = %d, longest length %d\n", typ, n, N, maxLen)

	// Fully commutative elements: R(x) commuting and xs fully commutative for
	// every s in R(x).
	fc := make([]bool, N)
	var fcCount uint32
	for x := uint32(0); x < N; x++ {
		ok := true
		r := Rmask[x]
		for i := 0; i < n && ok; i++ {
			if (r>>uint(i))&1 == 1 {
				for j := i + 1; j < n; j++ {
					if (r>>uint(j))&1 == 1 && M[i][j] > 2 {
						ok = false
						break
					}
				}
			}
		}
		for i := 0; i < n && ok; i++ {
			if (r>>uint(i))&1 == 1 && !fc[rmul[x][i]] {
				ok = false
			}
		}
		fc[x] = ok
		if ok {
			fcCount++
		}
	}

	// Reflections: one word per positive root (keys visited in lexicographic
	// order of the root vector, as std::map does).
	var reflWords [][]int
	{
		rootWord := map[[8]int][]int{}
		var q [][8]int
		for i := 0; i < n; i++ {
			var r [8]int
			r[i] = 1
			rootWord[r] = []int{i}
			q = append(q, r)
		}
		for h := 0; h < len(q); h++ {
			r := q[h]
			wd := rootWord[r]
			for i := 0; i < n; i++ {
				pair := 0
				for k := 0; k < n; k++ {
					pair += A[i][k] * r[k]
				}
				r2 := r
				r2[i] -= pair
				pos := true
				for k := 0; k < n; k++ {
					if r2[k] < 0 {
						pos = false
					}
				}
				if _, seen := rootWord[r2]; !pos || seen {
					continue
				}
				wd2 := make([]int, 0, len(wd)+2)
				wd2 = append(wd2, i)
				wd2 = append(wd2, wd...)
				wd2 = append(wd2, i)
				rootWord[r2] = wd2
				q = append(q, r2)
			}
		}
		keys := make([][8]int, 0, len(rootWord))
		for k := range rootWord {
			keys = append(keys, k)
		}
		sort.Slice(keys, func(a, b int) bool {
			for t := 0; t < 8; t++ {
				if keys[a][t] != keys[b][t] {
					return keys[a][t] < keys[b][t]
				}
			}
			return false
		})
		for _, k := range keys {
			reflWords = append(reflWords, rootWord[k])
		}
	}
	reflElt := make([]uint32, 0, len(reflWords))
	for _, wd := range reflWords {
		var e uint32
		for _, s := range wd {
			e = rmul[e][s]
		}
		reflElt = append(reflElt, e)
	}
	{
		distinct := map[uint32]bool{}
		for _, e := range reflElt {
			distinct[e] = true
		}
		if len(distinct) != len(reflWords) {
			fail(3, "reflections not distinct")
		}
	}

	// Bruhat order as bitsets: [e,w] = [e,v] union s[e,v] for w = sv, s in L(w).
	words := (int(N) + 63) / 64
	bruhat = make([][]uint64, N)
	for w := range bruhat {
		bruhat[w] = make([]uint64, words)
	}
	bruhat[0][0] = 1
	for w := uint32(1); w < N; w++ {
		s := bits.TrailingZeros8(Lmask[w])
		v := lmul[w][s]
		copy(bruhat[w], bruhat[v])
		for i := 0; i < words; i++ {
			b := bruhat[v][i]
			for b != 0 {
				t := bits.TrailingZeros64(b)
				b &= b - 1
				sx := lmul[uint32(i*64+t)][s]
				bruhat[w][sx>>6] |= 1 << (sx & 63)
			}
		}
	}
	var comparable uint64
	for w := uint32(0); w < N; w++ {
		for _, b := range bruhat[w] {
			comparable += uint64(bits.OnesCount64(b))
		}
	}
	// Deterministic cross-check of the bitsets against the lifting-property test.
	{
		state := uint64(88172645463325252)
		var checked uint64
		for t := 0; t < 200000; t++ {
			state ^= state << 13
			state ^= state >> 7
			state ^= state << 17
			x := uint32(state % uint64(N))
			state ^= state << 13
			state ^= state >> 7
			state ^= state << 17
			w := uint32(state % uint64(N))
			if leq(x, w) != leqIter(x, w) {
				fail(3, "Bruhat bitset mismatch")
			}
			checked++
		}
		fmt.Fprintf(os.Stderr, "Bruhat bitset and lifting comparisons match on %d sampled pairs\n", checked)
	}

	// Kazhdan-Lusztig recursion, by increasing length.
	ex = make([][]uint32, N)
	pid = make([][]uint32, N)
	mulist = make([][]zmu, N)
	internPoly([]int64{1})
	ex[0] = []uint32{0}
	pid[0] = []uint32{0}
	lenStart := make([]uint32, maxLen+2)
	for w := uint32(0); w < N; w++ {
		lenStart[slen[w]+1] = w + 1
	}
	for l := 1; l <= maxLen+1; l++ {
		if lenStart[l] == 0 {
			lenStart[l] = lenStart[l-1]
		}
	}
	extremalPairs := uint64(1)
	var cand []uint32
	var P, Q []int64
	resize := func(buf []int64, sz int) []int64 {
		if cap(buf) < sz {
			buf = make([]int64, sz)
		} else {
			buf = buf[:sz]
			for i := range buf {
				buf[i] = 0
			}
		}
		return buf
	}
	for l := 1; l <= maxLen; l++ {
		for w := lenStart[l]; w < lenStart[l+1]; w++ {
			s := bits.TrailingZeros8(Lmask[w])
			v := lmul[w][s]
			Lw, Rw := Lmask[w], Rmask[w]
			cand = cand[:0]
			for i := 0; i < words; i++ {
				bb := bruhat[w][i]
				for bb != 0 {
					t := bits.TrailingZeros64(bb)
					bb &= bb - 1
					x := uint32(i*64 + t)
					if Lmask[x]&Lw == Lw && Rmask[x]&Rw == Rw {
						cand = append(cand, x)
					}
				}
			}
			ids := make([]uint32, 0, len(cand))
			var mus []zmu
			for _, x := range cand {
				if x == w {
					ids = append(ids, 0)
					continue
				}
				d := int(slen[w]) - int(slen[x])
				P = resize(P, d/2+2)
				// Left recursion with s in L(w), w = s v; s in L(x) because x is extremal.
				p1, ok := lookup(lmul[x][s], v)
				if !ok {
					fail(3, "sx not below v")
				}
				for i := range p1 {
					P[i] += p1[i]
				}
				if p2, ok := lookup(x, v); ok {
					for i := range p2 {
						P[i+1] += p2[i]
					}
				}
				for _, zm := range mulist[v] {
					z := zm.z
					if (Lmask[z]>>uint(s))&1 == 0 || slen[z] < slen[x] {
						continue
					}
					p3, ok := lookup(x, z)
					if !ok {
						continue
					}
					shift := (int(slen[w]) - int(slen[z])) / 2
					for i := range p3 {
						P[i+shift] -= zm.mu * p3[i]
					}
				}
				guard(P)
				if P[0] != 1 {
					fail(3, "constant term != 1")
				}
				for _, c := range P {
					if c < 0 {
						fail(3, "negative coefficient")
					}
				}
				for i := (d-1)/2 + 1; i < len(P); i++ {
					if P[i] != 0 {
						fail(3, "degree bound violated")
					}
				}
				// Right recursion with t the largest right descent of w, w = v' t; t in R(x).
				{
					t := 7 - bits.LeadingZeros8(Rmask[w])
					v2 := rmul[w][t]
					Q = resize(Q, d/2+2)
					q1, ok := lookup(rmul[x][t], v2)
					if !ok {
						fail(3, "xt not below v'")
					}
					for i := range q1 {
						Q[i] += q1[i]
					}
					if q2, ok := lookup(x, v2); ok {
						for i := range q2 {
							Q[i+1] += q2[i]
						}
					}
					for _, zm := range mulist[v2] {
						z := zm.z
						if (Rmask[z]>>uint(t))&1 == 0 || slen[z] < slen[x] {
							continue
						}
						q3, ok := lookup(x, z)
						if !ok {
							continue
						}
						shift := (int(slen[w]) - int(slen[z])) / 2
						for i := range q3 {
							Q[i+shift] -= zm.mu * q3[i]
						}
					}
					if !eqPoly(P, Q) {
						fail(3, "left and right recursions return different polynomials")
					}
				}
				if d%2 == 1 {
					if mu := P[(d-1)/2]; mu != 0 {
						mus = append(mus, zmu{x, mu})
					}
				}
				ids = append(ids, internPoly(P))
			}
			// Bruhat covers that are not extremal have mu = 1.
			for ti := range reflElt {
				z := w
				for _, si := range reflWords[ti] {
					z = rmul[z][si]
				}
				if int(slen[z]) != int(slen[w])-1 {
					continue
				}
				if Lmask[z]&Lw == Lw && Rmask[z]&Rw == Rw {
					continue
				}
				mus = append(mus, zmu{z, 1})
			}
			sort.Slice(mus, func(a, b int) bool {
				if mus[a].z != mus[b].z {
					return mus[a].z < mus[b].z
				}
				return mus[a].mu < mus[b].mu
			})
			ex[w] = append([]uint32(nil), cand...)
			pid[w] = ids
			mulist[w] = mus
			extremalPairs += uint64(len(cand))
		}
	}
	var maxCoeff int64
	for _, p := range polys {
		for _, c := range p {
			if c > maxCoeff {
				maxCoeff = c
			}
		}
	}

	if *pairsFlag != "" {
		dumpPairs(*pairsFlag)
		return
	}

	// Statistics.
	muHist := map[int64]uint64{}
	var coverPairs, fcLowerNonzero, fcUpperNonzero uint64
	var maxMu, maxMuFcLower, maxMuFcUpper int64
	var big []bigEntry // mu >= 2
	for w := uint32(0); w < N; w++ {
		for _, zm := range mulist[w] {
			muHist[zm.mu]++
			if int(slen[w])-int(slen[zm.z]) == 1 {
				coverPairs++
			}
			if zm.mu > maxMu {
				maxMu = zm.mu
			}
			if zm.mu >= 2 {
				big = append(big, bigEntry{zm.z, w, zm.mu})
			}
			if fc[zm.z] {
				fcLowerNonzero++
				if zm.mu > maxMuFcLower {
					maxMuFcLower = zm.mu
				}
			}
			if fc[w] {
				fcUpperNonzero++
				if zm.mu > maxMuFcUpper {
					maxMuFcUpper = zm.mu
				}
			}
		}
	}
	wstr := map[uint32]string{}
	ws := func(e uint32) string {
		if s, ok := wstr[e]; ok {
			return s
		}
		s := wordStr(reducedWord(e))
		wstr[e] = s
		return s
	}
	sort.SliceStable(big, func(i, j int) bool {
		a, b := big[i], big[j]
		if a.mu != b.mu {
			return a.mu > b.mu
		}
		if slen[a.w] != slen[b.w] {
			return slen[a.w] < slen[b.w]
		}
		if wa, wb := ws(a.w), ws(b.w); wa != wb {
			return wa < wb
		}
		return ws(a.x) < ws(b.x)
	})
	var bigFcUpper uint64
	for _, t := range big {
		if fc[t.w] {
			bigFcUpper++
		}
	}

	pairJSON := func(x, w uint32, mu int64) string {
		p, _ := lookup(x, w)
		s := "    {\"mu\": " + strconv.FormatInt(mu, 10)
		s += ", \"x\": \"" + wordStr(reducedWord(x)) + "\", \"x_length\": " + strconv.Itoa(int(slen[x])) + ", \"x_fully_commutative\": " + boolStr(fc[x])
		s += ", \"x_L\": " + jsonInts(maskList(Lmask[x])) + ", \"x_R\": " + jsonInts(maskList(Rmask[x]))
		s += ", \"w\": \"" + wordStr(reducedWord(w)) + "\", \"w_length\": " + strconv.Itoa(int(slen[w])) + ", \"w_fully_commutative\": " + boolStr(fc[w])
		s += ", \"w_L\": " + jsonInts(maskList(Lmask[w])) + ", \"w_R\": " + jsonInts(maskList(Rmask[w]))
		s += ", \"gap\": " + strconv.Itoa(int(slen[w])-int(slen[x]))
		s += ", \"P_ascending\": " + jsonInts(p) + "}"
		return s
	}

	// Expected values for E6 (reviewer computation of 7 October 2026, recomputed here).
	e6 := typ == 'E' && n == 6
	if e6 {
		expect(N == 51840, "|W(E6)| = 51840")
		expect(maxLen == 36, "longest element of E6 has length 36")
		expect(fcCount == 662, "E6 has 662 fully commutative elements")
		expect(len(reflWords) == 36, "E6 has 36 positive roots")
		expect(maxMu == 10, "max mu over all pairs is 10")
		hist := func(k int64) uint64 { return muHist[k] }
		expect(hist(10) == 8, "8 pairs attain mu = 10")
		expect(hist(2) == 400 && hist(3) == 556 && hist(4) == 310 && hist(5) == 108 && hist(6) == 4, "histogram of mu >= 2 is 2:400 3:556 4:310 5:108 6:4 10:8")
		expect(hist(7) == 0 && hist(8) == 0 && hist(9) == 0, "no pair has mu in {7,8,9}")
		expect(len(big) == 1386, "1386 pairs have mu >= 2")
		expect(fcLowerNonzero == 6431 && maxMuFcLower == 1, "6431 pairs with FC lower endpoint and mu != 0, all with mu = 1")
		x, w := fromWord("01254012010"), fromWord("123012523412301251234012")
		expect(slen[x] == 11 && slen[w] == 24 && !fc[x] && !fc[w], "witness lengths 11 and 24, both non-FC")
		p, ok := lookup(x, w)
		expect(ok && eqPoly(p, []int64{1, 13, 59, 121, 125, 61, 10}), "P(01254012010, 123012523412301251234012) = 1+13q+59q^2+121q^3+125q^4+61q^5+10q^6")
		x5, w5 := fromWord("5343010"), fromWord("0125342312501234")
		expect(slen[x5] == 7 && slen[w5] == 16 && !fc[x5] && fc[w5], "FC-upper example: x of length 7 non-FC, w of length 16 FC")
		p5, ok5 := lookup(x5, w5)
		expect(ok5 && eqPoly(p5, []int64{1, 8, 22, 20, 5}), "P(5343010, 0125342312501234) = 1+8q+22q^2+20q^3+5q^4, so mu = 5")
		supp := map[int]bool{}
		for _, a := range reducedWord(w) {
			supp[a] = true
		}
		expect(len(supp) == 6, "witness w has full support")
	}

	var out bytes.Buffer
	fmt.Fprintf(&out, "{\n")
	fmt.Fprintf(&out, "  \"type\": \"%c%d\",\n", typ, n)
	labelling := "chain 0-1-...-(n-2), node n-1 attached to node 2; words act on the right, a string of labels is the product of the simple reflections in the order written"
	if typ == 'D' {
		labelling = "chain 0-...-(n-2), node n-1 attached to node n-3"
	} else if typ == 'A' {
		labelling = "chain 0-...-(n-1)"
	}
	fmt.Fprintf(&out, "  \"labelling\": \"%s\",\n", labelling)
	fmt.Fprintf(&out, "  \"group_order\": %d,\n  \"longest_length\": %d,\n  \"positive_roots\": %d,\n  \"fully_commutative_count\": %d,\n", N, maxLen, len(reflWords), fcCount)
	fmt.Fprintf(&out, "  \"bruhat_comparable_pairs\": %d,\n  \"extremal_pairs\": %d,\n  \"distinct_polynomials\": %d,\n  \"max_coefficient\": %d,\n", comparable, extremalPairs, len(polys), maxCoeff)
	fmt.Fprintf(&out, "  \"left_right_recursions_agree_on_every_extremal_pair\": true,\n")
	fmt.Fprintf(&out, "  \"max_mu\": %d,\n", maxMu)
	var nonzero uint64
	histKeys := make([]int64, 0, len(muHist))
	for k, v := range muHist {
		nonzero += v
		histKeys = append(histKeys, k)
	}
	sort.Slice(histKeys, func(i, j int) bool { return histKeys[i] < histKeys[j] })
	fmt.Fprintf(&out, "  \"nonzero_mu_pairs\": %d,\n  \"cover_pairs\": %d,\n", nonzero, coverPairs)
	fmt.Fprintf(&out, "  \"mu_histogram\": {")
	for i, k := range histKeys {
		sep := ", "
		if i == 0 {
			sep = ""
		}
		fmt.Fprintf(&out, "%s\"%d\": %d", sep, k, muHist[k])
	}
	fmt.Fprintf(&out, "},\n")
	fmt.Fprintf(&out, "  \"pairs_with_mu_at_least_2\": %d,\n", len(big))
	var bigFcLower uint64
	for _, b := range big {
		if fc[b.x] {
			bigFcLower++
		}
	}
	fmt.Fprintf(&out, "  \"pairs_with_mu_at_least_2_and_fully_commutative_lower_endpoint\": %d,\n", bigFcLower)
	fmt.Fprintf(&out, "  \"pairs_with_mu_at_least_2_and_fully_commutative_upper_endpoint\": %d,\n", bigFcUpper)
	fmt.Fprintf(&out, "  \"fully_commutative_lower_endpoint\": {\"nonzero_mu_pairs\": %d, \"max_mu\": %d},\n", fcLowerNonzero, maxMuFcLower)
	fmt.Fprintf(&out, "  \"fully_commutative_upper_endpoint\": {\"nonzero_mu_pairs\": %d, \"max_mu\": %d},\n", fcUpperNonzero, maxMuFcUpper)
	fmt.Fprintf(&out, "  \"pairs_attaining_max_mu\": [\n")
	{
		first := true
		for _, b := range big {
			if b.mu != maxMu {
				break
			}
			if !first {
				out.WriteString(",\n")
			}
			out.WriteString(pairJSON(b.x, b.w, b.mu))
			first = false
		}
	}
	fmt.Fprintf(&out, "\n  ],\n")
	fmt.Fprintf(&out, "  \"largest_mu_with_fully_commutative_upper_endpoint\": [\n")
	{
		first := true
		for _, b := range big {
			if b.mu != maxMuFcUpper || !fc[b.w] {
				continue
			}
			if !first {
				out.WriteString(",\n")
			}
			out.WriteString(pairJSON(b.x, b.w, b.mu))
			first = false
		}
	}
	fmt.Fprintf(&out, "\n  ],\n")
	if e6 {
		fmt.Fprintf(&out, "  \"expectations_checked\": \"E6 values of the manuscript\",\n")
	} else {
		fmt.Fprintf(&out, "  \"expectations_checked\": null,\n")
	}
	status := "passed"
	if failures != 0 {
		status = "FAILED"
	}
	fmt.Fprintf(&out, "  \"status\": \"%s\"\n}\n", status)
	if err := os.WriteFile(outfile, out.Bytes(), 0o644); err != nil {
		fail(2, "cannot open "+outfile)
	}
	fmt.Printf("%c%d: |W|=%d FC=%d extremal pairs=%d max mu=%d; pairs with mu>=2: %d; FC lower endpoint nonzero pairs %d with max mu %d; status %s\n",
		typ, n, N, fcCount, extremalPairs, maxMu, len(big), fcLowerNonzero, maxMuFcLower, status)
	if failures != 0 {
		os.Exit(1)
	}
}
