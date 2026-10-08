// Command terminal-structure certifies the structure of the noncommuting
// terminal elements of E6, E7 and E8 (the manuscript's remarks on the finite
// exceptional terminals): Gern's list plus two, the involution / layer /
// orthogonal-reflection structure of w_4, w_6, w_7, w_8, the right weak order
// chain w_4 < w_6 < w_7 < w_8, and the length-additive w_0(J) factorizations.
// Gern's Theorem 2.3.6 is re-derived for D5, D6 and D7 by exhaustive
// enumeration of signed permutations.
//
// It ports research/review_checks/terminal_structure.py. All arithmetic is
// exact (integer matrices; math/big.Rat for ranks).
//
// Imports: only the Go standard library. It imports no other package of this
// module, in particular none of the enumeration engines.
//
// Run (from the work directory): terminal-structure
// Writes results/terminal-structure-certificate.json and prints a one-line
// summary on stdout; exits non-zero if any expectation fails.
package main

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

type tableRow struct {
	n    int
	word string
	desc string
}

var table = []tableRow{ // (n, word, common descent set) as printed in Table 1
	{6, "1325213", "135"},
	{7, "1326213", "136"},
	{7, "13256213", "1356"},
	{7, "132543621324356", "1356"},
	{7, "1325436210321432543621324356", "1356"},
	{8, "1327213", "137"},
	{8, "13257213", "1357"},
	{8, "61327213", "1367"},
	{8, "132543721324357", "1357"},
	{8, "1325437210321432543721324357", "1357"},
	{8, "7534231270123456210321432" + "5437210321432543721324357", "1357"},
}

type layered struct {
	name   string
	n      int
	word   string
	layers [][]int
}

var layerOrder = []layered{ // layered palindromic words of the four maximal noncommuting terminals
	{"w4", 6, "1325213", [][]int{{1, 3}, {2}, {5}, {2}, {1, 3}}},
	{"w6", 7, "132543621324356", [][]int{{1, 3, 5, 6}, {2, 4}, {1, 3, 6}, {2, 4}, {1, 3, 5, 6}}},
	{"w7", 7, "1325436210321432543621324356",
		[][]int{{1, 3, 5, 6}, {2, 4}, {1, 3, 6}, {0, 2, 4}, {1, 3, 5, 6}, {0, 2, 4}, {1, 3, 6}, {2, 4}, {1, 3, 5, 6}}},
	{"w8", 8, "7534231270123456210321432" + "5437210321432543721324357",
		[][]int{{1, 3, 5, 7}, {2, 4}, {1, 3, 7}, {0, 2, 4}, {1, 3, 5, 7}, {0, 2, 4, 6}, {1, 3, 5, 7}, {2, 4},
			{1, 3, 5, 7}, {0, 2, 4, 6}, {1, 3, 5, 7}, {0, 2, 4}, {1, 3, 7}, {2, 4}, {1, 3, 5, 7}}},
}

var expectedNegated = map[string][]string{
	"w4": {"011101"}, "w6": {"0011001", "0111000", "0111111"},
	"w7": {"0011111", "0111101", "0121001", "1222111"}, "w8": {"12332212", "13432102"},
}

var expectedFactorization = map[string][3]string{
	"w7": {"31265234312", "653010", "42312645231"},
	"w8": {"31275234312", "753010", "276453423127563452341230127345231"},
}

func parse(word string) []int {
	out := make([]int, len(word))
	for i, c := range word {
		out[i] = int(c - '0')
	}
	return out
}

func wordStr(word []int) string {
	var sb strings.Builder
	for _, c := range word {
		sb.WriteString(strconv.Itoa(c))
	}
	return sb.String()
}

var failures []string

func expect(cond bool, what string) {
	if !cond {
		failures = append(failures, what)
		fmt.Fprintln(os.Stderr, "EXPECTATION FAILED:", what)
	}
}

func equalStrings(a, b []string) bool {
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

func equalInts(a, b []int) bool { return equalVec(a, b) }

func intsList(x []int) []int {
	if x == nil {
		return []int{}
	}
	return x
}

func strsList(x []string) []string {
	if x == nil {
		return []string{}
	}
	return x
}

func listRepr(x []int) string {
	parts := make([]string, len(x))
	for i, v := range x {
		parts[i] = strconv.Itoa(v)
	}
	return "[" + strings.Join(parts, ", ") + "]"
}

func partGernRows() obj {
	report := obj{}
	for _, n := range []int{4, 6, 8, 10} {
		word := gernWnWord(n)
		expect(equalInts(signedPerm(word, n), gernCor2219(n, n)), fmt.Sprintf("Lemma 2.3.4 word of w_%d equals Corollary 2.2.19", n))
		expect(len(word) == (3*n*n+2*n)/8, fmt.Sprintf("l(w_%d) = 3n^2/8 + n/4", n))
	}
	notFromD := map[int][]string{}
	for _, n := range []int{6, 7, 8} {
		rs := en(n)
		m := n - 1
		var rowOrder []string
		rows := map[string]mat{}
		for _, t := range table {
			if t.n == n {
				rowOrder = append(rowOrder, t.word)
				rows[t.word] = rs.wordMatrix(parse(t.word))
			}
		}
		matched := map[string]bool{}
		mapped := []any{}
		for _, g := range gernPredicted(m) {
			gword := append(gernWnWord(g.k), g.U...)
			mw := make([]int, len(gword))
			for i, s := range gword {
				mw[i] = gernToEn(s, n)
			}
			M := rs.wordMatrix(mw)
			var hits []string
			for _, w := range rowOrder {
				if equalMat(rows[w], M) {
					hits = append(hits, w)
				}
			}
			expect(len(hits) == 1, fmt.Sprintf("E%d: Gern %s matches exactly one table row", n, g.name))
			for _, h := range hits {
				matched[h] = true
			}
			var tr any
			if len(hits) > 0 {
				tr = hits[0]
			}
			mapped = append(mapped, obj{{"gern_element", g.name}, {"gern_word", wordStr(gword)},
				{"mapped_word", wordStr(mw)}, {"table_row", tr}, {"length", rs.length(M)}})
		}
		var unmatched []string
		for _, w := range rowOrder {
			if !matched[w] {
				unmatched = append(unmatched, w)
			}
		}
		sort.SliceStable(unmatched, func(a, b int) bool { return len(unmatched[a]) < len(unmatched[b]) })
		notFromD[n] = unmatched
		nodes := []int{1, n - 1}
		for x := 2; x < n-1; x++ {
			nodes = append(nodes, x)
		}
		sort.Ints(nodes)
		report = append(report, kv{fmt.Sprintf("E%d", n), obj{
			{"type_D_parabolic", fmt.Sprintf("D%d on nodes %s", m, listRepr(nodes))},
			{"gern_bad_elements_mapped", mapped},
			{"rows_not_from_type_D", strsList(unmatched)}}})
	}
	expect(len(notFromD[6]) == 0, "E6: every row is Gern's")
	expect(equalStrings(notFromD[7], []string{"1325436210321432543621324356"}), "E7: only w_7 is new")
	expect(equalStrings(notFromD[8], []string{"1325437210321432543721324357", table[len(table)-1].word}), "E8: only w_7 (in the E7 parabolic) and w_8 are not from type D")
	// The E8 row of length 28 is w_7 relabelled 6 -> 7 (E7 parabolic on {0,...,5,7}).
	rs8 := en(8)
	e7InE8 := strings.ReplaceAll(table[4].word, "6", "7")
	expect(equalMat(rs8.wordMatrix(parse(e7InE8)), rs8.wordMatrix(parse(table[9].word))), "E8 row of length 28 equals w_7 under 6 -> 7")
	report = append(report, kv{"E8_row_28_is_w7_in_parabolic_0_5_7", true})
	return report
}

func partDmEnumeration() obj {
	report := obj{}
	for _, m := range []int{5, 6, 7} {
		bad := enumerateBadDm(m)
		predicted := gernPredicted(m)
		predSet := map[string]bool{}
		for _, p := range predicted {
			predSet[keyInts(p.perm)] = true
		}
		badSet := map[string]bool{}
		for _, b := range bad {
			badSet[keyInts(b.w)] = true
		}
		equal := len(predSet) == len(badSet)
		for k := range badSet {
			if !predSet[k] {
				equal = false
			}
		}
		expect(equal, fmt.Sprintf("D%d: bad elements are exactly Gern's list", m))
		fact := 1
		for i := 2; i <= m; i++ {
			fact *= i
		}
		lengths := make([]int, len(bad))
		for i, b := range bad {
			lengths[i] = b.length
		}
		names := make([]string, len(predicted))
		for i, p := range predicted {
			names[i] = p.name
		}
		sort.Strings(names)
		report = append(report, kv{fmt.Sprintf("D%d", m), obj{
			{"elements", (1 << uint(m-1)) * fact}, {"noncommuting_terminals", len(bad)},
			{"lengths", intsList(lengths)}, {"gern_names", strsList(names)},
			{"equals_gern_theorem_2_3_6", equal}}})
	}
	return report
}

func containsInt(x []int, v int) bool {
	for _, a := range x {
		if a == v {
			return true
		}
	}
	return false
}

func sortedVecs(vs []vec) []vec {
	out := append([]vec(nil), vs...)
	sort.Slice(out, func(a, b int) bool { return cmpVec(out[a], out[b]) < 0 })
	return out
}

func rootStrs(rs []vec) []string {
	out := make([]string, len(rs))
	for i, r := range rs {
		out[i] = wordStr(r)
	}
	return out
}

func partStructure() obj {
	report := obj{}
	for _, lw := range layerOrder {
		name, n, wstr, layers := lw.name, lw.n, lw.word, lw.layers
		rs := en(n)
		W := rs.wordMatrix(parse(wstr))
		var full []int
		for _, layer := range layers {
			full = append(full, layer...)
		}
		expect(equalMat(rs.wordMatrix(full), W), name+": layered word gives the table element")
		expect(len(full) == rs.length(W) && rs.length(W) == len(wstr), name+": layered word is reduced")
		pal := true
		for i := range layers {
			if !equalInts(layers[i], layers[len(layers)-1-i]) {
				pal = false
			}
		}
		expect(pal, name+": layers are palindromic")
		colour := map[int]int{}
		for j := 0; j < n-1; j++ {
			colour[j] = j % 2
		}
		colour[n-1] = 1 // node n-1 is attached to node 2, so it has the colour of node 3
		alt := true
		classes := make([]map[int]bool, len(layers))
		for i, layer := range layers {
			classes[i] = map[int]bool{}
			for _, s := range layer {
				classes[i][colour[s]] = true
			}
			if len(classes[i]) != 1 {
				alt = false
			}
		}
		for i := 0; i+1 < len(classes); i++ {
			if len(classes[i]) == 1 && len(classes[i+1]) == 1 {
				for c := range classes[i] {
					if classes[i+1][c] {
						alt = false
					}
				}
			}
		}
		expect(alt, name+": layers alternate between the colour classes")
		indep := true
		for _, layer := range layers {
			for _, s := range layer {
				for _, t := range layer {
					if rs.adjSet[s][t] {
						indep = false
					}
				}
			}
		}
		expect(indep, name+": each layer is an independent set")
		L, R := rs.leftDescents(W), rs.rightDescents(W)
		tableDescents := map[string]bool{}
		for _, t := range table {
			if t.word == wstr {
				d := parse(t.desc)
				sort.Ints(d)
				tableDescents[keyInts(d)] = true
			}
		}
		expect(len(tableDescents) == 1 && tableDescents[keyInts(L)] && equalInts(L, R), name+": L = R = the common descent set of Table 1")
		sub := true
		for _, s := range layers[0] {
			if !containsInt(L, s) {
				sub = false
			}
		}
		expect(sub, name+": first layer is contained in L")
		expect(equalMat(rs.mult(W, W), rs.identity), name+": involution")
		expect(rs.terminal(W), name+": terminal")
		var negated, fixed []vec
		for _, r := range rs.positiveRoots() {
			img := rs.apply(W, r)
			neg := make(vec, len(r))
			for i, c := range r {
				neg[i] = -c
			}
			if equalVec(img, neg) {
				negated = append(negated, r)
			}
			if equalVec(img, r) {
				fixed = append(fixed, r)
			}
		}
		expect(equalStrings(rootStrs(negated), expectedNegated[name]), name+": negated positive roots")
		orth := true
		for i, a := range negated {
			for j, b := range negated {
				if i != j && rs.pairing(a, b) != 0 {
					orth = false
				}
			}
		}
		expect(orth, name+": negated roots mutually orthogonal")
		P := rs.identity
		for _, r := range negated {
			P = rs.mult(P, rs.reflection(r))
		}
		expect(equalMat(P, W), name+": product of the reflections in the negated roots")
		MI := make([][]int, n)
		for i := 0; i < n; i++ {
			MI[i] = make([]int, n)
			for j := 0; j < n; j++ {
				MI[i][j] = W[j][i]
				if i == j {
					MI[i][j]++
				}
			}
		}
		expect(n-rankOf(MI) == len(negated), name+": (-1)-eigenspace dimension equals the number of negated roots")
		// middle layer K and prefix P: W = P i(K) P^-1 with P(alpha_k), k in K, the negated roots
		half := layers[:len(layers)/2]
		K := layers[len(layers)/2]
		var halfWord []int
		for _, layer := range half {
			halfWord = append(halfWord, layer...)
		}
		Pm := rs.wordMatrix(halfWord)
		var images []vec
		for _, k := range K {
			e := make(vec, n)
			e[k] = 1
			images = append(images, rs.apply(Pm, e))
		}
		images = sortedVecs(images)
		sn := sortedVecs(negated)
		same := len(images) == len(sn)
		if same {
			for i := range images {
				if !equalVec(images[i], sn[i]) {
					same = false
				}
			}
		}
		expect(same, name+": images of the middle layer are the negated roots")
		entry := obj{{"type", fmt.Sprintf("E%d", n)}, {"word", wstr}, {"length", len(wstr)}, {"layers", layers},
			{"L_equals_R", intsList(L)}, {"first_layer", layers[0]},
			{"involution", true}, {"terminal", true}, {"negated_positive_roots", strsList(rootStrs(negated))},
			{"orthogonal_reflection_count", len(negated)}, {"fixed_positive_roots", len(fixed)}}
		if len(fixed) > 0 {
			rows := make([][]int, len(fixed))
			for i, r := range fixed {
				rows[i] = r
			}
			span := rankOf(rows)
			sub := obj{{"rank", span}, {"positive_roots", len(fixed)}}
			if name == "w8" {
				expect(span == 6 && len(fixed) == 30, "w8: fixed subsystem has rank 6 and 30 positive roots (type D6)")
				sub = append(sub, kv{"type", "D6"})
			}
			entry = append(entry, kv{"fixed_root_subsystem", sub})
		}
		report = append(report, kv{name, entry})
	}
	return report
}

func partChainAndFactorizations() obj {
	report := obj{}
	type pair struct{ a, b string }
	build := func(n int, words obj, pairs []pair, label string) obj {
		rs := en(n)
		mats := map[string]mat{}
		for _, e := range words {
			mats[e.K] = rs.wordMatrix(parse(e.V.(string)))
		}
		rel := []any{}
		for _, p := range pairs {
			pre, suf := rs.isPrefix(mats[p.a], mats[p.b]), rs.isSuffix(mats[p.a], mats[p.b])
			if label == "E8" {
				expect(pre && suf, fmt.Sprintf("E8: %s <= %s in right (and left) weak order", p.a, p.b))
			} else {
				expect(pre && suf, fmt.Sprintf("E7: %s <= %s in right weak order", p.a, p.b))
			}
			rel = append(rel, obj{{"lower", p.a}, {"upper", p.b}, {"prefix", pre}, {"suffix", suf}})
		}
		return obj{{"words", words}, {"relations", rel}}
	}
	// Chain in E8 and E7 (right weak order = prefix order; all four are involutions, so prefix = suffix).
	words8 := obj{{"w4", "1327213"}, {"w6", "132543721324357"}, {"w7", "1325437210321432543721324357"}, {"w8", table[len(table)-1].word}}
	words7 := obj{{"w4", "1326213"}, {"w6", "132543621324356"}, {"w7", "1325436210321432543621324356"}}
	e8 := build(8, words8, []pair{{"w4", "w6"}, {"w6", "w7"}, {"w7", "w8"}}, "E8")
	e7 := build(7, words7, []pair{{"w4", "w6"}, {"w6", "w7"}}, "E7")
	report = append(report, kv{"right_weak_order_chain", obj{{"E8", e8}, {"E7", e7}}})

	// Longest-parabolic inner factors and the w_0(J) factorizations.
	factors := obj{}
	maxLen := map[string]int{}
	for _, lw := range layerOrder {
		name, n, wstr := lw.name, lw.n, lw.word
		rs := en(n)
		W := rs.wordMatrix(parse(wstr))
		suff := rs.suffixes(W)
		best := 0
		var bestJ []int
		for _, V := range suff.order {
			Lv := rs.leftDescents(V)
			if len(Lv) > 0 {
				l := rs.length(rs.longestElement(Lv))
				if l > best {
					best, bestJ = l, Lv
				}
			}
		}
		maxLen[name] = best
		entry := obj{{"type", fmt.Sprintf("E%d", n)}, {"right_weak_interval_size", len(suff.order)},
			{"max_inner_longest_parabolic_length", best}, {"achieved_by_J", intsList(bestJ)}}
		I := rs.leftDescents(W)
		Jset := map[int]bool{0: true}
		for _, s := range I {
			Jset[s] = true
		}
		var J []int
		for s := range Jset {
			J = append(J, s)
		}
		sort.Ints(J)
		w0J := rs.longestElement(J)
		var hit mat
		for _, V := range suff.order {
			ld := rs.leftDescents(V)
			all := true
			for _, s := range J {
				if !containsInt(ld, s) {
					all = false
				}
			}
			if all {
				hit = V
				break
			}
		}
		if hit != nil {
			u := suff.u[keyMat(hit)]
			v := rs.mult(rs.inverse(w0J), hit)
			lu, lw0, lv := rs.length(u), rs.length(w0J), rs.length(v)
			additive := lu+lw0+lv == rs.length(W)
			expect(equalMat(rs.mult(rs.mult(u, w0J), v), W) && additive, name+": w = u w_0(J) v length-additively")
			var jEdges [][]int
			for _, e := range rs.edges {
				if Jset[e[0]] && Jset[e[1]] {
					jEdges = append(jEdges, []int{e[0], e[1]})
				}
			}
			uw, w0w, vw := wordStr(rs.reducedWord(u)), wordStr(rs.reducedWord(w0J)), wordStr(rs.reducedWord(v))
			entry = append(entry, kv{"w0_factorization", obj{{"J", J}, {"l_w0_J", lw0}, {"u", uw}, {"w0_J", w0w}, {"v", vw},
				{"lengths", []int{lu, lw0, lv}}, {"length_additive", additive}, {"J_edges", jEdges}}})
			if exp, ok := expectedFactorization[name]; ok {
				expect(equalMat(rs.wordMatrix(parse(exp[0])), u) && equalMat(rs.wordMatrix(parse(exp[1])), w0J) && equalMat(rs.wordMatrix(parse(exp[2])), v),
					name+": factorization matches the manuscript's words")
				expect(lw0 == 6 && len(jEdges) == 1, name+": J = I + {0} is of type A2 x A1^3 with l(w_0(J)) = 6")
			}
		}
		factors = append(factors, kv{name, entry})
	}
	expect(maxLen["w4"] == 3, "w4: longest inner parabolic factor has length 3")
	expect(maxLen["w6"] == 5, "w6: longest inner parabolic factor has length 5")
	expect(maxLen["w7"] == 6, "w7: longest inner parabolic factor has length 6")
	expect(maxLen["w8"] == 6, "w8: longest inner parabolic factor has length 6")
	report = append(report, kv{"inner_longest_parabolic_factors", factors})
	return report
}

func get(o obj, k string) any {
	for _, e := range o {
		if e.K == k {
			return e.V
		}
	}
	panic("missing key " + k)
}

func main() {
	report := obj{
		{"conventions", "E_n: chain 0-1-...-(n-2), node n-1 attached to node 2; a string of labels is the product of the simple reflections in the order written; columns of the matrix of w are w(alpha_j). Gern D_m labels: 1,2 attached to 3, chain 3-...-m; label map 1->1, 2->n-1, j->j-1."},
		{"gern_plus_two", partGernRows()},
		{"type_D_enumeration", partDmEnumeration()},
		{"structure", partStructure()},
	}
	report = append(report, partChainAndFactorizations()...)
	status := "passed"
	if len(failures) > 0 {
		status = "FAILED"
	}
	report = append(report, kv{"status", status}, kv{"failed_expectations", strsList(failures)})

	outPath := filepath.Join("results", "terminal-structure-certificate.json")
	if err := os.MkdirAll("results", 0o755); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if err := os.WriteFile(outPath, []byte(toJSON(report, 1)+"\n"), 0o644); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}

	var dRanks []string
	for _, e := range get(report, "type_D_enumeration").(obj) {
		dRanks = append(dRanks, e.K)
	}
	sort.Strings(dRanks)
	var chain8 []string
	for _, r := range get(get(get(report, "right_weak_order_chain").(obj), "E8").(obj), "relations").([]any) {
		ro := r.(obj)
		chain8 = append(chain8, get(ro, "lower").(string)+"<"+get(ro, "upper").(string))
	}
	inner := obj{}
	for _, e := range get(report, "inner_longest_parabolic_factors").(obj) {
		inner = append(inner, kv{e.K, get(e.V.(obj), "max_inner_longest_parabolic_length")})
	}
	fmt.Println(toJSON(obj{{"status", status}, {"type_D_ranks_enumerated", dRanks}, {"chain_E8", chain8},
		{"inner_factor_lengths", inner}}, -1))
	if len(failures) > 0 {
		os.Exit(1)
	}
}
