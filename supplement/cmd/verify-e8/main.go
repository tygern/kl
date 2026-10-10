// Command verify-e8 is the independent integral-matrix audit of the two E8
// coset certificates (the E7 / 240-coset and the D7 / 2160-coset terminal
// classifications) together with the complete FC catalogue.
//
// It ports research/en_independent/verify_e8.py. It checks root inversion
// lengths, word equality across the D7 and E7 decompositions, complete
// FC-catalogue closure, all eligible bottoms and the exact parabolic/coset
// Poincare products, and writes results/e8-independent-audit.json (the
// summary is also printed on stdout). Unlike the original it hashes only its
// JSON inputs: the Go supplement carries no Python or C++ source, and the
// archive's MANIFEST.json hashes every Go source file.
//
// Imports: only the Go standard library. It imports no other package of this
// module, in particular none of the enumeration engines (terminals-flat,
// terminals-recursive, e8-d7, fc-catalogue): all arithmetic (integer column
// matrices, root orbit, inversion-count lengths) is local to this package.
//
// Inputs are read relative to the current working directory (the work
// directory of the runner). Run: verify-e8
package main

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math/big"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

const n = 8

// edges is the E8 Dynkin diagram: the chain 0-1-...-6 plus the node 7 on 2.
var edges = [][2]int{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {2, 7}}

var adj [n][]int

var identity elem

// vec is an integer vector of length 8 (a root or a matrix column).
type vec [n]int64

// elem is the integer matrix of a Weyl group element, stored as its columns.
type elem [n]vec

func init() {
	for _, e := range edges {
		adj[e[0]] = append(adj[e[0]], e[1])
		adj[e[1]] = append(adj[e[1]], e[0])
	}
	for j := 0; j < n; j++ {
		identity[j][j] = 1
	}
}

// failf reports a failed assertion and exits non-zero.
func failf(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "verify-e8: assertion failed: "+format+"\n", args...)
	os.Exit(1)
}

func assert(cond bool, format string, args ...any) {
	if !cond {
		failf(format, args...)
	}
}

// addCheck, mulCheck guard int64 accumulation against overflow.
func addCheck(a, b int64) int64 {
	c := a + b
	if (c > a) != (b > 0) {
		failf("int64 overflow in addition")
	}
	return c
}

func mulCheck(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}
	c := a * b
	if c/b != a {
		failf("int64 overflow in multiplication")
	}
	return c
}

func neg(v vec) vec {
	var r vec
	for i := range v {
		r[i] = -v[i]
	}
	return r
}

func plus(u, v vec) vec {
	var r vec
	for i := range u {
		r[i] = u[i] + v[i]
	}
	return r
}

// right multiplies the matrix a by the simple reflection s on the right.
func right(a elem, s int) elem {
	cols := a
	cols[s] = neg(a[s])
	for _, t := range adj[s] {
		cols[t] = plus(a[t], a[s])
	}
	return cols
}

func element(word []int) elem {
	a := identity
	for _, s := range word {
		assert(s >= 0 && s < n, "generator %d out of range", s)
		a = right(a, s)
	}
	return a
}

func reversed(word []int) []int {
	r := make([]int, len(word))
	for i, s := range word {
		r[len(word)-1-i] = s
	}
	return r
}

// positive reports whether the vector is a positive root; it asserts that the
// vector is nonzero with all coordinates of one sign.
func positive(a vec) bool {
	anyNonzero, anyPos, anyNeg := false, false, false
	for _, v := range a {
		if v != 0 {
			anyNonzero = true
		}
		if v > 0 {
			anyPos = true
		}
		if v < 0 {
			anyNeg = true
		}
	}
	assert(anyNonzero && !(anyPos && anyNeg), "column is not a root: %v", a)
	return anyPos
}

// desc returns the right descent set as a bit mask (bit s set iff column s of
// a is a negative root). mask(desc(a)) of the original is the value itself.
func desc(a elem) uint {
	var m uint
	for s := 0; s < n; s++ {
		if !positive(a[s]) {
			m |= 1 << uint(s)
		}
	}
	return m
}

func independent(ds uint) bool {
	for _, e := range edges {
		if ds&(1<<uint(e[0])) != 0 && ds&(1<<uint(e[1])) != 0 {
			return false
		}
	}
	return true
}

// weak tests that no suffix ts with s,t adjacent exists, without words.
func weak(a elem) bool {
	d := desc(a)
	for s := 0; s < n; s++ {
		if d&(1<<uint(s)) == 0 {
			continue
		}
		for _, t := range adj[s] {
			if !positive(plus(a[s], a[t])) {
				return false
			}
		}
	}
	return true
}

// image computes sum_s root[s]*a[s] (coordinatewise).
func image(a elem, root vec) vec {
	var r vec
	for k := 0; k < n; k++ {
		var sum int64
		for s := 0; s < n; s++ {
			sum = addCheck(sum, mulCheck(root[s], a[s][k]))
		}
		r[k] = sum
	}
	return r
}

func roots() map[vec]bool {
	found := map[vec]bool{}
	var todo []vec
	for j := 0; j < n; j++ {
		found[identity[j]] = true
		todo = append(todo, identity[j])
	}
	for i := 0; i < len(todo); i++ {
		a := todo[i]
		for s := 0; s < n; s++ {
			b := a
			sum := int64(0)
			for _, t := range adj[s] {
				sum += a[t]
			}
			b[s] = -a[s] + sum
			if !found[b] {
				found[b] = true
				todo = append(todo, b)
			}
		}
	}
	assert(len(found) == 240, "root count %d, expected 240", len(found))
	for r := range found {
		assert(positive(r) || positive(neg(r)), "root %v has mixed signs", r)
	}
	return found
}

func minBit(m uint) int {
	if m == 0 {
		failf("min of empty descent set")
	}
	for s := 0; s < n; s++ {
		if m&(1<<uint(s)) != 0 {
			return s
		}
	}
	return -1
}

// leq decides x <= w in the Bruhat order by the descent-subword recursion.
func leq(x elem, lx int, w elem, lw int) bool {
	for x != w {
		if lx >= lw {
			return false
		}
		s := minBit(desc(w))
		if desc(x)&(1<<uint(s)) != 0 {
			x, lx = right(x, s), lx-1
		}
		w, lw = right(w, s), lw-1
	}
	return true
}

var volatileKeys = []string{"seconds", "elapsed_seconds"}

func isVolatile(k string) bool {
	for _, v := range volatileKeys {
		if k == v {
			return true
		}
	}
	return false
}

// writeCanonical serialises a decoded JSON value exactly like Python's
// json.dumps(sort_keys=True, separators=(',', ':')) with ensure_ascii, after
// removing the volatile keys at every depth. Non-integer numbers are
// rejected (none occur outside the volatile keys).
func writeCanonical(buf *bytes.Buffer, v any) {
	switch x := v.(type) {
	case nil:
		buf.WriteString("null")
	case bool:
		if x {
			buf.WriteString("true")
		} else {
			buf.WriteString("false")
		}
	case json.Number:
		bi, ok := new(big.Int).SetString(x.String(), 10)
		assert(ok, "non-integer number %s in hashed JSON", x.String())
		buf.WriteString(bi.String())
	case string:
		writeString(buf, x)
	case []any:
		buf.WriteByte('[')
		for i, e := range x {
			if i > 0 {
				buf.WriteByte(',')
			}
			writeCanonical(buf, e)
		}
		buf.WriteByte(']')
	case map[string]any:
		keys := make([]string, 0, len(x))
		for k := range x {
			if !isVolatile(k) {
				keys = append(keys, k)
			}
		}
		sort.Strings(keys)
		buf.WriteByte('{')
		for i, k := range keys {
			if i > 0 {
				buf.WriteByte(',')
			}
			writeString(buf, k)
			buf.WriteByte(':')
			writeCanonical(buf, x[k])
		}
		buf.WriteByte('}')
	default:
		failf("unexpected JSON value type %T", v)
	}
}

func writeString(buf *bytes.Buffer, s string) {
	buf.WriteByte('"')
	for _, r := range s {
		switch {
		case r == '"':
			buf.WriteString(`\"`)
		case r == '\\':
			buf.WriteString(`\\`)
		case r == '\n':
			buf.WriteString(`\n`)
		case r == '\r':
			buf.WriteString(`\r`)
		case r == '\t':
			buf.WriteString(`\t`)
		case r == '\b':
			buf.WriteString(`\b`)
		case r == '\f':
			buf.WriteString(`\f`)
		case r < 0x20 || r == 0x7f || (r >= 0x80 && r < 0x10000):
			fmt.Fprintf(buf, `\u%04x`, r)
		case r >= 0x10000:
			r -= 0x10000
			fmt.Fprintf(buf, `\u%04x\u%04x`, 0xd800+(r>>10), 0xdc00+(r&0x3ff))
		default:
			buf.WriteRune(r)
		}
	}
	buf.WriteByte('"')
}

func decodeGeneric(data []byte) any {
	dec := json.NewDecoder(bytes.NewReader(data))
	dec.UseNumber()
	var v any
	if err := dec.Decode(&v); err != nil {
		failf("invalid JSON: %v", err)
	}
	return v
}

func contentHash(path string) string {
	data, err := os.ReadFile(path)
	if err != nil {
		failf("cannot read %s: %v", path, err)
	}
	if strings.HasSuffix(path, ".json") {
		var buf bytes.Buffer
		writeCanonical(&buf, decodeGeneric(data))
		sum := sha256.Sum256(buf.Bytes())
		return hex.EncodeToString(sum[:])
	}
	sum := sha256.Sum256(data)
	return hex.EncodeToString(sum[:])
}

func convolution(a, b []int64) []int64 {
	c := make([]int64, len(a)+len(b)-1)
	for i, x := range a {
		for j, y := range b {
			c[i+j] = addCheck(c[i+j], mulCheck(x, y))
		}
	}
	return c
}

func poincare(exponents []int) []int64 {
	a := []int64{1}
	for _, m := range exponents {
		ones := make([]int64, m+1)
		for i := range ones {
			ones[i] = 1
		}
		a = convolution(a, ones)
	}
	return a
}

func equalInts(a, b []int64) bool {
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

type bottom struct {
	Word   []int `json:"word"`
	Length int   `json:"length"`
}

type badRow struct {
	Word    []int    `json:"word"`
	Length  int      `json:"length"`
	Rmask   uint     `json:"Rmask"`
	Lmask   uint     `json:"Lmask"`
	Bottoms []bottom `json:"bottoms"`
}

type fcRow struct {
	Word   []int `json:"word"`
	Length int   `json:"length"`
	Rmask  uint  `json:"Rmask"`
	Lmask  uint  `json:"Lmask"`
}

type d7Cert struct {
	Complete           bool     `json:"complete"`
	ParabolicOrder     int64    `json:"parabolic_order"`
	Cosets             int64    `json:"cosets"`
	ParabolicHistogram []int64  `json:"parabolic_histogram"`
	CosetHistogram     []int64  `json:"coset_histogram"`
	FullRightTerminals int64    `json:"full_right_terminals"`
	CommutingTerminals int64    `json:"commuting_terminals"`
	Bad                []badRow `json:"bad"`
}

type e7Cert struct {
	Complete           bool     `json:"complete"`
	ParabolicOrder     int64    `json:"parabolic_order"`
	CosetCount         int64    `json:"coset_count"`
	RightTerminalCount int64    `json:"right_terminal_count"`
	CommutingTerminals int64    `json:"commuting_terminals"`
	Bad                []badRow `json:"bad"`
}

type catalogueCert struct {
	Complete      bool    `json:"complete"`
	AscentsTested int64   `json:"ascents_tested"`
	Elements      []fcRow `json:"elements"`
}

// loadTyped decodes the file into out after checking that every required key
// is present (a missing key is a KeyError in the original).
func loadTyped(path string, out any, required []string) {
	data, err := os.ReadFile(path)
	if err != nil {
		failf("cannot read %s: %v", path, err)
	}
	var keys map[string]json.RawMessage
	if err := json.Unmarshal(data, &keys); err != nil {
		failf("%s: invalid JSON: %v", path, err)
	}
	for _, k := range required {
		_, ok := keys[k]
		assert(ok, "%s: missing key %q", path, k)
	}
	if err := json.Unmarshal(data, out); err != nil {
		failf("%s: %v", path, err)
	}
}

type record struct {
	length int
	rd, ld uint
	elig   elem
	eligL  int
}

type fcItem struct {
	a    elem
	info fcRow
}

type kv struct{ k, v string }

// sourceHashes marshals as a JSON object preserving insertion order.
type sourceHashes []kv

func (h sourceHashes) MarshalJSON() ([]byte, error) {
	var buf bytes.Buffer
	buf.WriteByte('{')
	for i, e := range h {
		if i > 0 {
			buf.WriteByte(',')
		}
		k, _ := json.Marshal(e.k)
		v, _ := json.Marshal(e.v)
		buf.Write(k)
		buf.WriteByte(':')
		buf.Write(v)
	}
	buf.WriteByte('}')
	return buf.Bytes(), nil
}

type summary struct {
	Status                           string       `json:"status"`
	CompleteE8TerminalClassification bool         `json:"complete_E8_terminal_classification"`
	DistinctParabolicMethods         []string     `json:"distinct_parabolic_methods"`
	MatrixWordsAndLengthsAgree       bool         `json:"matrix_words_and_inversion_lengths_agree"`
	TerminalCount                    int          `json:"terminal_count"`
	CommutingTerminals               int          `json:"commuting_terminals"`
	BadTerminalLengths               []int        `json:"bad_terminal_lengths"`
	EligibleRanks                    []int        `json:"eligible_ranks"`
	CompleteFCCatalogue              int          `json:"complete_FC_catalogue"`
	FCAscentsVerified                int          `json:"FC_ascents_verified"`
	PoincareProductsMatch            bool         `json:"parabolic_and_coset_poincare_products_match"`
	SourceSHA256                     sourceHashes `json:"source_sha256"`
	ProgramSources                   string       `json:"program_sources"`
	Hashing                          string       `json:"hashing"`
	Priority                         string       `json:"priority"`
}

func main() {
	jsonPaths := []string{
		"research/en_independent/e8-d7-terminals.json",
		"research/en_e8/e8-terminals.json",
		"research/en_independent/e8-fc.json",
	}
	// The original also hashed its own source and the three generator
	// sources (verify_e8.py, e8_d7_cosets.cpp, fc_catalogue.cpp,
	// parabolic_terminals.cpp) byte for byte. The Go supplement ships no
	// Python or C++ source, and every Go source file is hashed by the
	// archive's MANIFEST.json, so this audit records only the JSON inputs
	// and names the Go programs whose outputs it audits.
	programSources := "cmd/verify-e8 (this audit), cmd/e8-d7, cmd/fc-catalogue, cmd/terminals-flat and internal/parabolic; their Go source files are hashed by the supplement's MANIFEST.json, not here."
	var d7 d7Cert
	var e7 e7Cert
	var catalogue catalogueCert
	loadTyped(jsonPaths[0], &d7, []string{"complete", "parabolic_order", "cosets", "parabolic_histogram", "coset_histogram", "full_right_terminals", "commuting_terminals", "bad"})
	loadTyped(jsonPaths[1], &e7, []string{"complete", "parabolic_order", "coset_count", "right_terminal_count", "commuting_terminals", "bad"})
	loadTyped(jsonPaths[2], &catalogue, []string{"complete", "ascents_tested", "elements"})
	assert(d7.Complete && e7.Complete && catalogue.Complete, "a certificate is not marked complete")
	assert(mulCheck(d7.ParabolicOrder, d7.Cosets) == 696729600, "D7 parabolic order times cosets")
	assert(mulCheck(e7.ParabolicOrder, e7.CosetCount) == 696729600, "E7 parabolic order times cosets")
	assert(equalInts(d7.ParabolicHistogram, poincare([]int{1, 3, 5, 7, 9, 11, 6})), "D7 parabolic histogram")
	assert(equalInts(convolution(d7.ParabolicHistogram, d7.CosetHistogram), poincare([]int{1, 7, 11, 13, 17, 19, 23, 29})), "D7 Poincare product")
	assert(d7.FullRightTerminals == 2160 && e7.RightTerminalCount == 2160, "right terminal counts")
	assert(d7.CommutingTerminals == 58 && e7.CommutingTerminals == 58, "commuting terminal counts")
	indep := 0
	for bits := uint(0); bits < 1<<n; bits++ {
		if independent(bits) {
			indep++
		}
	}
	assert(indep == 58, "independent subsets: %d", indep)

	allRoots := roots()
	var positiveRoots []vec
	for r := range allRoots {
		if positive(r) {
			positiveRoots = append(positiveRoots, r)
		}
	}
	length := func(a elem) int {
		c := 0
		for _, r := range positiveRoots {
			if !positive(image(a, r)) {
				c++
			}
		}
		return c
	}

	fc := map[elem]fcRow{}
	var fcOrder []fcItem
	for _, row := range catalogue.Elements {
		a := element(row.Word)
		_, dup := fc[a]
		assert(!dup && independent(desc(a)), "FC element %v duplicate or with non-independent descents", row.Word)
		assert(len(row.Word) == row.Length, "FC word length mismatch %v", row.Word)
		assert(desc(a) == row.Rmask, "FC Rmask mismatch %v", row.Word)
		assert(desc(element(reversed(row.Word))) == row.Lmask, "FC Lmask mismatch %v", row.Word)
		fc[a] = row
		fcOrder = append(fcOrder, fcItem{a, row})
	}
	assert(len(fc) == 10846, "FC catalogue size %d", len(fc))

	// By the recursive FC criterion and a minimal-missing-element argument,
	// the finite catalogue contains exactly the FC elements.
	ascents := 0
	for _, it := range fcOrder {
		a, row := it.a, it.info
		ds := desc(a)
		for s := 0; s < n; s++ {
			if ds&(1<<uint(s)) == 0 {
				continue
			}
			lower, ok := fc[right(a, s)]
			assert(ok, "descent of FC element not in catalogue")
			assert(lower.Length == row.Length-1, "descent length drop")
		}
		for s := 0; s < n; s++ {
			if ds&(1<<uint(s)) != 0 {
				continue
			}
			ascents++
			b := right(a, s)
			db := desc(b)
			valid := independent(db)
			if valid {
				for t := 0; t < n; t++ {
					if db&(1<<uint(t)) == 0 {
						continue
					}
					if _, ok := fc[right(b, t)]; !ok {
						valid = false
						break
					}
				}
			}
			_, inFC := fc[b]
			assert(valid == inFC, "recursive FC test and catalogue membership differ")
		}
	}
	assert(int64(ascents) == catalogue.AscentsTested, "ascents %d vs %d", ascents, catalogue.AscentsTested)

	var canonical [2]map[elem]record
	for idx, bad := range [][]badRow{d7.Bad, e7.Bad} {
		records := map[elem]record{}
		for _, row := range bad {
			a := element(row.Word)
			ai := element(reversed(row.Word))
			_, dup := records[a]
			_, inFC := fc[a]
			assert(!dup && weak(a) && weak(ai) && !inFC, "bad element %v duplicate, not weak, or in FC", row.Word)
			assert(length(a) == row.Length && row.Length == len(row.Word), "bad element length %v", row.Word)
			images := map[vec]bool{}
			for r := range allRoots {
				images[image(a, r)] = true
			}
			ok := len(images) == len(allRoots)
			for r := range allRoots {
				if !images[r] {
					ok = false
				}
			}
			assert(ok, "bad element %v does not permute the roots", row.Word)
			rd, ld := desc(a), desc(ai)
			assert(rd == row.Rmask && ld == row.Lmask, "bad element masks %v", row.Word)
			type pair struct {
				x elem
				l int
			}
			var eligible []pair
			for _, it := range fcOrder {
				info := it.info
				if info.Rmask&row.Rmask != row.Rmask || info.Lmask&row.Lmask != row.Lmask {
					continue
				}
				if leq(it.a, info.Length, a, row.Length) {
					eligible = append(eligible, pair{it.a, info.Length})
				}
			}
			supplied := map[pair]bool{}
			for _, b := range row.Bottoms {
				supplied[pair{element(b.Word), b.Length}] = true
			}
			elSet := map[pair]bool{}
			for _, e := range eligible {
				elSet[e] = true
			}
			same := len(elSet) == len(supplied)
			for e := range elSet {
				if !supplied[e] {
					same = false
				}
			}
			assert(same, "eligible bottoms differ from supplied bottoms for %v", row.Word)
			assert(len(eligible) == 1, "bad element %v has %d eligible bottoms", row.Word, len(eligible))
			records[a] = record{row.Length, rd, ld, eligible[0].x, eligible[0].l}
		}
		canonical[idx] = records
	}
	same := len(canonical[0]) == len(canonical[1])
	for a, r0 := range canonical[0] {
		if r1, ok := canonical[1][a]; !ok || r0 != r1 {
			same = false
		}
	}
	assert(same, "D7 and E7 terminal records differ")

	var lengths, ranks []int
	for _, r := range canonical[0] {
		lengths = append(lengths, r.length)
		ranks = append(ranks, r.length-r.eligL)
	}
	sort.Ints(lengths)
	sort.Ints(ranks)

	var hashes sourceHashes
	for _, p := range jsonPaths {
		hashes = append(hashes, kv{filepath.ToSlash(p), contentHash(p)})
	}
	vk := append([]string{}, volatileKeys...)
	sort.Strings(vk)
	out := summary{
		Status:                           "passed",
		CompleteE8TerminalClassification: true,
		DistinctParabolicMethods:         []string{"E7, 240 cosets", "D7, 2160 cosets"},
		MatrixWordsAndLengthsAgree:       true,
		TerminalCount:                    64,
		CommutingTerminals:               58,
		BadTerminalLengths:               lengths,
		EligibleRanks:                    ranks,
		CompleteFCCatalogue:              len(fc),
		FCAscentsVerified:                ascents,
		PoincareProductsMatch:            true,
		SourceSHA256:                     hashes,
		ProgramSources:                   programSources,
		Hashing:                          "JSON inputs are hashed after removing the volatile keys " + strings.Join(vk, ", ") + " and re-serialising with sorted keys.",
		Priority:                         "Not established by this computational audit.",
	}
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(out); err != nil {
		failf("encoding summary: %v", err)
	}
	if err := os.MkdirAll("results", 0o755); err != nil {
		failf("%v", err)
	}
	if err := os.WriteFile("results/e8-independent-audit.json", buf.Bytes(), 0o644); err != nil {
		failf("%v", err)
	}
	os.Stdout.Write(buf.Bytes())
}
