// Command verify-exceptional compares the independent finite terminal
// certificates of the exceptional groups E6 and E7 and checks them against
// known invariants.
//
// Four snapshots are compared: the root-index enumerations of E6 and E7
// (research/broad_exceptional/e6_bad.json and e7_bad.json, both written by the
// enumerate-bad command) and the integer-matrix enumerations
// (results/e6-independent-certificate.json written by matrix-search -rank 6,
// results/e7-independent-certificate.json written by matrix-search -rank 7).
// The two programs of each pair use different representations and share no
// code. This program only reads the four JSON snapshots; it never runs any
// enumeration engine and does not need exact_e.py (historical provenance).
//
// Ports research/verify_exceptional.py (audit mode; the --recompute mode is
// replaced by the runner executing enumerate-bad and matrix-search). Reads the
// four snapshots relative to the current directory and writes
// results/exceptional-audit.json (also printed on stdout).
//
// Imports: the Go standard library only. It imports no internal package and
// none of the other engines.
package main

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
)

func fail(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "verify-exceptional: FAILED: "+format+"\n", args...)
	os.Exit(1)
}

func require(cond bool, format string, args ...any) {
	if !cond {
		fail(format, args...)
	}
}

// ---- JSON access helpers (a missing key or wrong type is a failure, as the
// Python KeyError/TypeError would be) ----

func readFile(path string) []byte {
	data, err := os.ReadFile(path)
	if err != nil {
		fail("cannot read %s: %v", path, err)
	}
	return data
}

func readJSON(path string) map[string]any {
	dec := json.NewDecoder(bytes.NewReader(readFile(path)))
	dec.UseNumber()
	var v any
	if err := dec.Decode(&v); err != nil {
		fail("cannot parse %s: %v", path, err)
	}
	if _, err := dec.Token(); err != io.EOF {
		fail("%s: trailing data after JSON value", path)
	}
	m, ok := v.(map[string]any)
	if !ok {
		fail("%s: top level is not an object", path)
	}
	return m
}

func field(m map[string]any, key string) any {
	v, ok := m[key]
	if !ok {
		fail("missing key %q", key)
	}
	return v
}

func asInt(v any, what string) int64 {
	n, ok := v.(json.Number)
	if !ok {
		fail("%s: not a number", what)
	}
	i, err := strconv.ParseInt(n.String(), 10, 64)
	if err != nil {
		fail("%s: not an integer: %v", what, err)
	}
	return i
}

func asList(v any, what string) []any {
	l, ok := v.([]any)
	if !ok {
		fail("%s: not a list", what)
	}
	return l
}

func asMap(v any, what string) map[string]any {
	m, ok := v.(map[string]any)
	if !ok {
		fail("%s: not an object", what)
	}
	return m
}

func intField(m map[string]any, key string) int64 { return asInt(field(m, key), key) }

func intList(v any, what string) []int64 {
	l := asList(v, what)
	out := make([]int64, len(l))
	for i, x := range l {
		out[i] = asInt(x, what)
	}
	return out
}

// ---- canonical rows (tuples compared lexicographically, like Python) ----

type bottom struct {
	word   []int64
	length int64
	rank   int64
}

type row struct {
	word    []int64
	length  int64
	bottoms []bottom
}

func cmpInts(a, b []int64) int {
	for i := 0; i < len(a) && i < len(b); i++ {
		if a[i] != b[i] {
			if a[i] < b[i] {
				return -1
			}
			return 1
		}
	}
	switch {
	case len(a) < len(b):
		return -1
	case len(a) > len(b):
		return 1
	}
	return 0
}

func cmpInt(a, b int64) int {
	if a < b {
		return -1
	}
	if a > b {
		return 1
	}
	return 0
}

func cmpBottom(a, b bottom) int {
	if c := cmpInts(a.word, b.word); c != 0 {
		return c
	}
	if c := cmpInt(a.length, b.length); c != 0 {
		return c
	}
	return cmpInt(a.rank, b.rank)
}

func cmpRow(a, b row) int {
	if c := cmpInts(a.word, b.word); c != 0 {
		return c
	}
	if c := cmpInt(a.length, b.length); c != 0 {
		return c
	}
	for i := 0; i < len(a.bottoms) && i < len(b.bottoms); i++ {
		if c := cmpBottom(a.bottoms[i], b.bottoms[i]); c != 0 {
			return c
		}
	}
	return cmpInt(int64(len(a.bottoms)), int64(len(b.bottoms)))
}

func canonical(rows any, what, bottomsKey, rankKey string) []row {
	list := asList(rows, what)
	out := make([]row, len(list))
	for i, r := range list {
		m := asMap(r, what)
		bs := asList(field(m, bottomsKey), what+"."+bottomsKey)
		bots := make([]bottom, len(bs))
		for j, b := range bs {
			bm := asMap(b, what+"."+bottomsKey)
			bots[j] = bottom{intList(field(bm, "word"), "word"), intField(bm, "length"), intField(bm, rankKey)}
		}
		out[i] = row{intList(field(m, "word"), "word"), intField(m, "length"), bots}
	}
	sort.SliceStable(out, func(i, j int) bool { return cmpRow(out[i], out[j]) < 0 })
	return out
}

// Rows of enumerate-bad: word, length, Rmask, Lmask, bottoms[word, length, rank].
func canonicalRootIndex(rows any, what string) []row {
	return canonical(rows, what, "bottoms", "rank")
}

// Rows of the E6 matrix program: word, length, candidates[word, length, interval_rank].
func canonicalE6Matrix(rows any, what string) []row {
	return canonical(rows, what, "candidates", "interval_rank")
}

// Rows of the E7 matrix program: word, length, bottoms[word, length, rank].
func canonicalE7Matrix(rows any, what string) []row {
	return canonicalRootIndex(rows, what)
}

func equalRows(a, b []row) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if cmpRow(a[i], b[i]) != 0 {
			return false
		}
	}
	return true
}

func equalInts(a, b []int64) bool { return cmpInts(a, b) == 0 }

// polynomial multiplies the factors 1 + q + ... + q^exponent.
func polynomial(exponents []int) []int64 {
	result := []int64{1}
	for _, exponent := range exponents {
		next := make([]int64, len(result)+exponent)
		for i, c := range result {
			for j := 0; j <= exponent; j++ {
				next[i+j] += c
				require(next[i+j] >= 0, "overflow in polynomial")
			}
		}
		result = next
	}
	return result
}

func sum(xs []int64) int64 {
	var s int64
	for _, x := range xs {
		s += x
		require(s >= 0, "overflow in sum")
	}
	return s
}

func signed(word []int64) []int64 {
	mapping := map[int64]int{1: 0, 6: 1, 2: 2, 3: 3, 4: 4, 5: 5}
	values := []int64{1, 2, 3, 4, 5, 6}
	for _, label := range word {
		s, ok := mapping[label]
		require(ok, "unknown generator label %d in signed permutation check", label)
		if s == 0 {
			values[0], values[1] = -values[1], -values[0]
		} else {
			values[s-1], values[s] = values[s], values[s-1]
		}
	}
	return values
}

func fileSHA(path string) string {
	h := sha256.Sum256(readFile(path))
	return hex.EncodeToString(h[:])
}

type generators struct {
	E6Bad    string `json:"research/broad_exceptional/e6_bad.json"`
	E7Bad    string `json:"research/broad_exceptional/e7_bad.json"`
	E6Matrix string `json:"results/e6-independent-certificate.json"`
	E7Matrix string `json:"results/e7-independent-certificate.json"`
}

// The original also hashed the four generator sources; those programs are not
// part of this module's payload, so only the four data snapshots are hashed.
type fileHashes struct {
	E6Matrix string `json:"results/e6-independent-certificate.json"`
	E7Matrix string `json:"results/e7-independent-certificate.json"`
	E6Bad    string `json:"research/broad_exceptional/e6_bad.json"`
	E7Bad    string `json:"research/broad_exceptional/e7_bad.json"`
}

type report struct {
	Status                         string     `json:"status"`
	IndependentTerminalTablesAgree bool       `json:"independent_terminal_tables_agree"`
	PoincareHistogramsMatch        bool       `json:"poincare_histograms_match"`
	D6SignedPairMatches            bool       `json:"d6_signed_pair_matches"`
	SnapshotGenerators             generators `json:"snapshot_generators"`
	Files                          fileHashes `json:"files"`
	Priority                       string     `json:"priority"`
}

func main() {
	if len(os.Args) > 1 {
		fail("takes no arguments")
	}
	const (
		pE6Bad    = "research/broad_exceptional/e6_bad.json"
		pE7Bad    = "research/broad_exceptional/e7_bad.json"
		pE6Matrix = "results/e6-independent-certificate.json"
		pE7Matrix = "results/e7-independent-certificate.json"
	)
	a6 := readJSON(pE6Bad)
	b6 := readJSON(pE6Matrix)
	a7 := readJSON(pE7Bad)
	b7 := readJSON(pE7Matrix)

	require(equalRows(canonicalRootIndex(field(a6, "bad"), "e6 bad"), canonicalE6Matrix(field(b6, "non_fc_weak_bad_terminals"), "e6 matrix bad")),
		"E6 root-index and matrix terminal tables differ")
	require(equalRows(canonicalRootIndex(field(a7, "bad"), "e7 bad"), canonicalE7Matrix(field(b7, "bad"), "e7 matrix bad")),
		"E7 root-index and matrix terminal tables differ")

	// Poincare polynomials of E6 and E7 from their exponents, element counts,
	// fully commutative counts (Stembridge 1998) and commuting terminal counts.
	e6 := polynomial([]int{1, 4, 5, 7, 8, 11})
	e7 := polynomial([]int{1, 5, 7, 9, 11, 13, 17})
	require(sum(e6) == 51840 && sum(e7) == 2903040, "Poincare polynomial totals wrong")

	a6dist := intList(field(a6, "length_distribution"), "a6 length_distribution")
	b6distMap := asMap(field(b6, "length_distribution"), "b6 length_distribution")
	b6dist := make([]int64, len(e6))
	for i := range e6 {
		b6dist[i] = asInt(field(b6distMap, strconv.Itoa(i)), "b6 length_distribution entry")
	}
	require(equalInts(a6dist, b6dist) && equalInts(b6dist, e6), "E6 length distributions differ from the Poincare polynomial")
	a7dist := intList(field(a7, "length_distribution"), "a7 length_distribution")
	b7dist := intList(field(b7, "length_distribution"), "b7 length_distribution")
	require(equalInts(a7dist, b7dist) && equalInts(b7dist, e7), "E7 length distributions differ from the Poincare polynomial")

	require(intField(a6, "order") == intField(b6, "element_count") && intField(b6, "element_count") == 51840, "E6 order mismatch")
	require(intField(a7, "order") == intField(b7, "order") && intField(b7, "order") == 2903040, "E7 order mismatch")
	require(intField(a6, "fc_count") == intField(b6, "fc_count") && intField(b6, "fc_count") == 662, "E6 fc_count mismatch")
	require(intField(a7, "fc_count") == intField(b7, "fc_count") && intField(b7, "fc_count") == 2670, "E7 fc_count mismatch")
	require(intField(a6, "commuting_weak") == intField(b6, "commuting_terminals") && intField(b6, "commuting_terminals") == 22, "E6 commuting count mismatch")
	require(intField(a7, "commuting_weak") == intField(b7, "commuting_terminals") && intField(b7, "commuting_terminals") == 36, "E7 commuting count mismatch")

	// Star-operation closure checks: the matrix programs count both sides.
	require(intField(b6, "fc_star_checks") == 2*intField(a6, "fc_right_star_checks") && 2*intField(a6, "fc_right_star_checks") == 3620, "E6 star-check count mismatch")
	require(intField(b7, "fc_star_checks") == 2*intField(a7, "fc_right_star_checks") && 2*intField(a7, "fc_right_star_checks") == 16776, "E7 star-check count mismatch")

	// Independent signed permutation check of the odd terminal, without a KL
	// evaluator or imports from either exceptional enumeration.
	var odd map[string]any
	for _, r := range asList(field(b7, "bad"), "b7 bad") {
		m := asMap(r, "b7 bad row")
		if intField(m, "length") == 15 {
			odd = m
			break
		}
	}
	require(odd != nil, "no E7 bad terminal of length 15")
	require(equalInts(signed(intList(field(odd, "word"), "odd word")), []int64{-1, -6, 3, -4, 5, -2}), "signed permutation of the odd terminal wrong")
	oddBottoms := asList(field(odd, "bottoms"), "odd bottoms")
	require(len(oddBottoms) > 0, "odd terminal has no bottoms")
	bottomWord := intList(field(asMap(oddBottoms[0], "odd bottom"), "word"), "odd bottom word")
	require(equalInts(signed(bottomWord), []int64{-1, -2, 4, 3, 6, 5}), "signed permutation of the odd bottom wrong")

	rep := report{
		Status:                         "passed",
		IndependentTerminalTablesAgree: true,
		PoincareHistogramsMatch:        true,
		D6SignedPairMatches:            true,
		SnapshotGenerators: generators{
			E6Bad:    "enumerate-bad -rank 6",
			E7Bad:    "enumerate-bad -rank 7",
			E6Matrix: "matrix-search -rank 6",
			E7Matrix: "matrix-search -rank 7",
		},
		Files: fileHashes{
			E6Matrix: fileSHA(pE6Matrix),
			E7Matrix: fileSHA(pE7Matrix),
			E6Bad:    fileSHA(pE6Bad),
			E7Bad:    fileSHA(pE7Bad),
		},
		Priority: "Not established by this computational audit.",
	}
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(rep); err != nil {
		fail("encoding report: %v", err)
	}
	if err := os.WriteFile("results/exceptional-audit.json", buf.Bytes(), 0o644); err != nil {
		fail("writing results/exceptional-audit.json: %v", err)
	}
	os.Stdout.Write(buf.Bytes())
}
