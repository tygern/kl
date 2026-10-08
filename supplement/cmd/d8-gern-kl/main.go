// Command d8-gern-kl computes Kazhdan-Lusztig polynomials on the Bruhat lower
// ideal of Gern's bad element w_n in W(D_n), n even, for n = 6 and n = 8, in
// two independent models of the group:
//
//	model sp:  signed permutations, Gern's conventions (thesis Example 1.1.6,
//	           Proposition 2.2.4): s_1 acts on the right by (w_1,w_2) -> (-w_2,-w_1),
//	           s_i (i >= 2) swaps positions i-1 and i;
//	model geo: the geometric representation on the root lattice (Cartan matrix
//	           of D_n with Gern's labels: nodes 1 and 2 attached to 3, chain 3-4-...-n).
//
// The engine stores only extremal pairs (x,y) with R(y) subset R(x) and
// L(y) subset L(x) and uses the right-descent recurrence (Kazhdan-Lusztig 1979,
// (2.2.c) transported to the right).  At the top element it recomputes
// P_{x,w} for every x in the ideal by the direct recurrence for every right
// descent of w, and all values must agree with the stored ones.
//
// Certified: D6 (P_{x_6,w_6} = 1+6q+11q^2+6q^3+q^4+q^5) and D8
// (P_{x_8,w_8} = P_{e,w_8} = 1+12q+59q^2+154q^3+233q^4+221q^5+147q^6+70q^7+20q^8+2q^9,
// mu(x_8,w_8) = 0).  Gern's element w_n is taken from Corollary 2.2.19
// (signed permutation) and checked against the bracket word of Lemma 2.3.4;
// x_n = s_1 s_2 s_4 s_6 ... s_n.
//
// Ports research/review_checks/d8_gern_kl.cpp.  Standard library only; imports
// no other package of this module (self-contained review program).
//
// Concurrency: the elements of one length layer are distributed over
// goroutines; each goroutine writes only its own result slots and reads only
// tables of lower layers.  Polynomials are interned sequentially in layer
// order after the join, so all results are independent of scheduling and of
// the thread count.
//
// Usage: d8-gern-kl [-out results/d8-gern-certificate.json] [-threads 0] [-models sp|geo|both]
//
// Hidden debugging flag (differential testing only, no certificate output):
// d8-gern-kl -pairs FILE [-pairsgeo] reads lines "n wword xword" (n = rank
// of D_n, 4 <= n <= 8; words as strings of Gern's labels 1..n, "e" for the
// identity), builds the Bruhat lower ideal of w in the signed-permutation
// model, runs the same extremal-pair engine, and prints one line
// "n wword xword P_sp P_geo" with P_{x,w} in ascending powers ("[]" when x is
// not below w); P_geo is computed by the geometric model when -pairsgeo is
// given and is "-" otherwise.  With -pairs absent the program's behaviour and
// output are unchanged.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"math/bits"
	"os"
	"runtime"
	"strconv"
	"strings"
)

var failures = 0

func expect(cond bool, what string) {
	if !cond {
		fmt.Fprintf(os.Stderr, "EXPECTATION FAILED: %s\n", what)
		failures++
	}
}

type result struct {
	ideal, lx, lw, interval int
	rankIdeal, rankInterval []int64
	Px, Pe                  poly
	mu                      int64
	mismatches              int
	polys                   int
	maxc                    int64
	canonTop, canonX        string
}

func runModel(model string, n, nth int) result {
	var topword, xword []int // Gern labels 1..n
	br := func(j, i int) []int {
		var w []int
		if i == 0 {
			for k := j; k >= 1; k-- {
				w = append(w, k)
			}
		} else {
			for k := j; k >= i; k-- {
				w = append(w, k)
			}
		}
		return w
	}
	{
		k := n/2 - 2
		for j := 2; j <= n; j += 2 {
			topword = append(topword, br(j, 0)...)
		}
		for i := 0; i <= k; i++ {
			topword = append(topword, br(n-k+i, n-2*k+2*i)...)
		}
	}
	xword = []int{1, 2}
	for j := 4; j <= n; j += 2 {
		xword = append(xword, j)
	}
	var M *Model
	if model == "sp" {
		var top perm
		for i := 1; i <= n; i++ {
			var v int
			switch {
			case i == 1:
				if (n/2)%2 != 0 {
					v = -1
				} else {
					v = 1
				}
			case i%2 != 0:
				v = i
			default:
				v = -(n + 2 - i)
			}
			top[i-1] = int8(v) // Corollary 2.2.19
		}
		M = spBuild(n, top, xword)
		g := spGroup{n}
		w := g.identity()
		for _, s := range topword {
			w = g.rmul(w, s-1)
		}
		if !g.equal(w, top) {
			fatal("Lemma 2.3.4 word and Corollary 2.2.19 disagree")
		}
		if g.length(top) != (3*n*n+2*n)/8 {
			fatal("length of w_n is not 3n^2/8 + n/4")
		}
	} else {
		M = geoBuild(n, topword, xword)
	}
	if M.idX < 0 {
		fatal("x_n is not in the lower ideal of w_n")
	}
	if len(topword) != int(M.length[M.idTop]) {
		fatal("bracket word not reduced")
	}
	canonicalWords(M)
	kl := newKL(M)
	kl.run(nth)
	top := M.idTop
	rd := M.rdes[top]
	var R result
	R.ideal = M.N
	R.lx = int(M.length[M.idX])
	R.lw = int(M.length[top])
	R.canonTop = M.canon[top]
	R.canonX = M.canon[M.idX]
	R.rankIdeal = make([]int64, R.lw+1)
	R.rankInterval = make([]int64, R.lw-R.lx+1)
	for x := 0; x < M.N; x++ {
		var ref poly
		first := true
		for s := 0; s < n; s++ {
			if (rd>>uint(s))&1 != 0 {
				P := kl.compute(x, top, s, int(M.rmul[top][s]))
				if first {
					ref = P
					first = false
				} else if P != ref {
					R.mismatches++
				}
			}
		}
		var extP poly
		if e := kl.lookup(x, top); e != nil {
			extP = *e
		}
		if extP != ref {
			R.mismatches++
		}
		R.rankIdeal[M.length[x]]++
	}
	for z := 0; z < M.N; z++ {
		if kl.lookup(M.idX, z) != nil {
			R.rankInterval[int(M.length[z])-R.lx]++
			R.interval++
		}
	}
	px := kl.lookup(M.idX, top)
	if px == nil {
		fatal("P(x,w) missing from the table")
	}
	R.Px = *px
	s0 := bits.TrailingZeros8(rd)
	R.Pe = kl.compute(M.idE, top, s0, int(M.rmul[top][s0]))
	cod := R.lw - R.lx
	if cod%2 == 0 {
		R.mu = 0
	} else {
		R.mu = R.Px[(cod-1)/2]
	}
	R.polys = len(kl.polys)
	for j := range kl.polys {
		p := &kl.polys[j]
		if p[0] != 1 {
			fatal("constant term != 1")
		}
		for i := 0; i < maxDeg; i++ {
			if p[i] < 0 {
				fatal("negative coefficient")
			}
			if p[i] > R.maxc {
				R.maxc = p[i]
			}
		}
	}
	fmt.Fprintf(os.Stderr, "D%d %s: ideal %d, l(w)=%d, l(x)=%d, interval %d, P(x,w)=%s, mismatches %d\n",
		n, model, R.ideal, R.lw, R.lx, R.interval, polyJSON(&R.Px), R.mismatches)
	return R
}

// parseGernWord reads a word of Gern labels 1..n ("e" for the identity) for
// the -pairs debugging flag and returns zero-based generator indices.
func parseGernWord(w string, n int) []int {
	if w == "e" {
		return nil
	}
	out := make([]int, 0, len(w))
	for _, c := range []byte(w) {
		if c < '1' || int(c-'0') > n {
			fatal("bad Gern label in word " + w)
		}
		out = append(out, int(c-'1'))
	}
	return out
}

// pairModel caches one built ideal (and its KL table) per top element.
type pairModel struct {
	M  *Model
	kl *KL
}

func buildPairModel(M *Model) *pairModel {
	kl := newKL(M)
	kl.run(1)
	return &pairModel{M, kl}
}

// walk follows a reduced word of x through the right multiplication table
// from the identity; -1 means that x is outside the ideal (every prefix of a
// reduced word of x lies below x, so a prefix outside the ideal shows that x
// is not below the top).
func (pm *pairModel) walk(word []int) int {
	id := pm.M.idE
	for _, s := range word {
		id = int(pm.M.rmul[id][s])
		if id < 0 {
			return -1
		}
	}
	return id
}

func (pm *pairModel) poly(x int) string {
	if x < 0 {
		return "[]"
	}
	p := pm.kl.lookup(x, pm.M.idTop)
	if p == nil {
		return "[]"
	}
	return polyJSON(p)
}

// runPairs implements the hidden -pairs flag (see the file header).  Each
// line "n wword xword" is answered by "n wword xword P_sp P_geo".
func runPairs(path string, withGeo bool) {
	f, err := os.Open(path)
	if err != nil {
		fatal("cannot open pairs file " + path)
	}
	defer f.Close()
	sp := map[string]*pairModel{}
	geo := map[string]*pairModel{}
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		fields := strings.Fields(sc.Text())
		if len(fields) == 0 || strings.HasPrefix(fields[0], "#") {
			continue
		}
		if len(fields) < 3 {
			fatal("pairs line needs rank and two words")
		}
		n, err := strconv.Atoi(fields[0])
		if err != nil || n < 4 || n > 8 {
			fatal("rank must be between 4 and 8")
		}
		key := fields[0] + " " + fields[1]
		topword := parseGernWord(fields[1], n)
		xword := parseGernWord(fields[2], n)
		g := spGroup{n}
		// Canonical reduced word of x (strip the smallest right descent).
		x := g.identity()
		for _, s := range xword {
			x = g.rmul(x, s)
		}
		var xred []int
		for u := x; g.length(u) > 0; {
			for s := 0; s < n; s++ {
				if v := g.rmul(u, s); g.length(v) < g.length(u) {
					xred = append(xred, s)
					u = v
					break
				}
			}
		}
		for i, j := 0, len(xred)-1; i < j; i, j = i+1, j-1 {
			xred[i], xred[j] = xred[j], xred[i]
		}
		pm, ok := sp[key]
		if !ok {
			top := g.identity()
			for _, s := range topword {
				top = g.rmul(top, s)
			}
			pm = buildPairModel(spBuild(n, top, nil))
			sp[key] = pm
		}
		pGeo := "-"
		if withGeo {
			gm, ok := geo[key]
			if !ok {
				gern := make([]int, len(topword))
				for i, s := range topword {
					gern[i] = s + 1
				}
				gm = buildPairModel(geoBuild(n, gern, nil))
				geo[key] = gm
			}
			pGeo = gm.poly(gm.walk(xred))
		}
		fmt.Printf("%s %s %s %s %s\n", fields[0], fields[1], fields[2], pm.poly(pm.walk(xred)), pGeo)
	}
	if err := sc.Err(); err != nil {
		fatal("reading pairs file: " + err.Error())
	}
}

func eqVec(a, b []int64) bool {
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

func vecJSON(v []int64) string {
	parts := make([]string, len(v))
	for i, c := range v {
		parts[i] = fmt.Sprint(c)
	}
	return "[" + strings.Join(parts, ",") + "]"
}

func main() {
	out := flag.String("out", "results/d8-gern-certificate.json", "output certificate path")
	threads := flag.Int("threads", 0, "worker goroutines (0 = min(8, number of CPUs))")
	models := flag.String("models", "sp", "group models to run: sp, geo or both")
	pairs := flag.String("pairs", "", "debugging: file of \"n wword xword\" lines; print P_{x,w} per line and write no certificate")
	pairsGeo := flag.Bool("pairsgeo", false, "debugging: with -pairs, also compute P_{x,w} in the geometric model")
	flag.Parse()
	if flag.NArg() != 0 {
		fmt.Fprintln(os.Stderr, "usage: d8-gern-kl [-out file] [-threads n] [-models sp|geo|both]")
		os.Exit(2)
	}
	if *pairs != "" {
		runPairs(*pairs, *pairsGeo)
		return
	}
	if *models != "sp" && *models != "geo" && *models != "both" {
		fmt.Fprintln(os.Stderr, "-models must be sp, geo or both")
		os.Exit(2)
	}
	nth := *threads
	if nth <= 0 {
		nth = runtime.NumCPU()
		if nth > 8 {
			nth = 8
		}
		if nth < 1 {
			nth = 1
		}
	}
	which := *models

	var sb strings.Builder
	sb.WriteString("{\n  \"group_conventions\": \"W(D_n), Gern's labels 1..n: s_1 acts on the right by (w_1,w_2)->(-w_2,-w_1), s_i (i>=2) swaps positions i-1,i; nodes 1 and 2 attached to 3, chain 3-...-n\",\n")
	sb.WriteString("  \"models\": [\"signed permutations\", \"geometric representation on the root lattice\"],\n  \"ranks\": {\n")
	expectedD6 := []int64{1, 6, 11, 6, 1, 1}
	expectedD8 := []int64{1, 12, 59, 154, 233, 221, 147, 70, 20, 2}
	for _, n := range []int{6, 8} {
		first := "sp"
		if which == "geo" {
			first = "geo"
		}
		a := runModel(first, n, nth)
		b := a
		if which == "both" {
			b = runModel("geo", n, nth)
		}
		expect(a.Px == b.Px && a.Pe == b.Pe && a.ideal == b.ideal && a.interval == b.interval &&
			eqVec(a.rankIdeal, b.rankIdeal) && eqVec(a.rankInterval, b.rankInterval), "the two models agree")
		expect(a.mismatches == 0 && b.mismatches == 0, "direct recurrence for every right descent agrees with the stored table")
		px := a.Px[:deg(&a.Px)+1]
		pe := a.Pe[:deg(&a.Pe)+1]
		if n == 6 {
			expect(eqVec(px, expectedD6), "P(x_6,w_6) = 1+6q+11q^2+6q^3+q^4+q^5")
			expect(a.lw == 15 && a.lx == 4 && a.ideal == 3184 && a.interval == 1676 && a.mu == 1, "D6: l(w)=15, l(x)=4, ideal 3184, interval 1676, mu = 1")
			expect(eqVec(a.rankInterval, []int64{1, 12, 55, 140, 248, 339, 360, 287, 162, 59, 12, 1}), "D6 interval rank vector")
		} else {
			expect(eqVec(px, expectedD8) && eqVec(pe, expectedD8), "P(x_8,w_8) = P(e,w_8) = 1+12q+59q^2+154q^3+233q^4+221q^5+147q^6+70q^7+20q^8+2q^9")
			expect(a.lw == 26 && a.lx == 5 && a.ideal == 265760 && a.interval == 163724 && a.mu == 0, "D8: l(w)=26, l(x)=5, ideal 265,760, interval 163,724, mu = 0")
		}
		fmt.Fprintf(&sb, "    \"D%d\": {\n", n)
		fmt.Fprintf(&sb, "      \"w_canonical_word\": \"%s\", \"x_canonical_word\": \"%s\",\n", a.canonTop, a.canonX)
		fmt.Fprintf(&sb, "      \"w_length\": %d, \"x_length\": %d, \"lower_ideal_size\": %d, \"interval_size\": %d,\n", a.lw, a.lx, a.ideal, a.interval)
		fmt.Fprintf(&sb, "      \"lower_ideal_rank_vector\": %s,\n", vecJSON(a.rankIdeal))
		fmt.Fprintf(&sb, "      \"interval_rank_vector\": %s,\n", vecJSON(a.rankInterval))
		fmt.Fprintf(&sb, "      \"P_x_w_ascending\": %s,\n      \"P_e_w_ascending\": %s,\n      \"degree\": %d,\n      \"mu_x_w\": %d,\n", polyJSON(&a.Px), polyJSON(&a.Pe), deg(&a.Px), a.mu)
		sep := ""
		if n == 6 {
			sep = ","
		}
		fmt.Fprintf(&sb, "      \"distinct_polynomials\": %d, \"max_coefficient\": %d,\n      \"models_agree\": true, \"direct_recurrence_mismatches\": %d\n    }%s\n", a.polys, a.maxc, a.mismatches+b.mismatches, sep)
	}
	status := "passed"
	if failures != 0 {
		status = "FAILED"
	}
	fmt.Fprintf(&sb, "  },\n  \"x_n\": \"s_1 s_2 s_4 s_6 ... s_n\",\n  \"w_n\": \"Gern Corollary 2.2.19, equal to the bracket word of Lemma 2.3.4\",\n  \"status\": \"%s\"\n}\n", status)
	if err := os.WriteFile(*out, []byte(sb.String()), 0o644); err != nil {
		fmt.Fprintln(os.Stderr, "cannot open output:", err)
		os.Exit(2)
	}
	fmt.Printf("D6 and D8 Gern-element Kazhdan-Lusztig certificate: %s\n", status)
	if failures != 0 {
		os.Exit(1)
	}
}
