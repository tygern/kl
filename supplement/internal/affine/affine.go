// Package affine is the shared affine-E8 arithmetic used by cmd/affine-proof
// and cmd/affine-fc-covers.
//
// Its row integer matrices, 240 finite E8 roots and inversion formulas port
// the library part of research/en_affine_referee/verify_families.py. The family
// audit certifies reduced translation prefixes and explicit descent braids
// without an FC catalogue. cmd/affine-fc-covers separately uses the catalogue
// for the supplementary closure and complete base-cover checks.
//
// Standard library only. No KL computation and no affine group enumeration.
// Any failed assertion of the original is a panic carrying an
// AssertionError; commands convert it into a non-zero exit via Main.
package affine

import (
	"bytes"
	"encoding/json"
	"fmt"
	"os"
	"sort"
)

// N is the number of generators (affine E8 with nine nodes).
const N = 9

// Vec is an integer row vector; Mat is a row matrix.
type Vec [N]int64
type Mat [N][N]int64

// limit guards against silent int64 overflow; real entries stay far smaller.
const limit = int64(1) << 28

// AssertionError is the panic value used for every failed assertion.
type AssertionError struct{ Msg string }

func (e AssertionError) Error() string { return "assertion failed: " + e.Msg }

// Assert panics with an AssertionError when cond is false.
func Assert(cond bool, format string, args ...any) {
	if !cond {
		panic(AssertionError{fmt.Sprintf(format, args...)})
	}
}

// Main runs f and exits with status 1 and a clear message on any panic.
func Main(f func()) {
	defer func() {
		if r := recover(); r != nil {
			fmt.Fprintf(os.Stderr, "FAILED: %v\n", r)
			os.Exit(1)
		}
	}()
	f()
}

// Edges is the E9 = affine E8 diagram.
var Edges = func() [][2]int {
	e := [][2]int{}
	for i := 0; i < 7; i++ {
		e = append(e, [2]int{i, i + 1})
	}
	return append(e, [2]int{2, 8})
}()

// Adj holds the sorted neighbours of each node.
var Adj = func() [N][]int {
	var adj [N][]int
	for _, e := range Edges {
		adj[e[0]] = append(adj[e[0]], e[1])
		adj[e[1]] = append(adj[e[1]], e[0])
	}
	for i := range adj {
		sort.Ints(adj[i])
	}
	return adj
}()

// E is the identity matrix.
var E = func() Mat {
	var m Mat
	for i := 0; i < N; i++ {
		m[i][i] = 1
	}
	return m
}()

var (
	Delta = Vec{2, 4, 6, 5, 4, 3, 2, 1, 3}
	Beta  = Vec{1, 2, 3, 3, 2, 2, 1, 1, 2}
	Gamma = Vec{1, 1, 1, 1, 1, 0, 0, 0, 0}
)

func abs(x int64) int64 {
	if x < 0 {
		return -x
	}
	return x
}

// Add returns a + k*b.
func Add(a, b Vec, k int64) Vec {
	var out Vec
	for i := range out {
		out[i] = a[i] + k*b[i]
		Assert(abs(out[i]) <= limit, "coefficient bound exceeded")
	}
	return out
}

// Dot is the Euclidean dot product.
func Dot(a, b Vec) int64 {
	var s int64
	for i := range a {
		Assert(abs(a[i]) <= limit && abs(b[i]) <= limit, "coefficient bound exceeded")
		s += a[i] * b[i]
	}
	return s
}

// Pair is the symmetric bilinear form of the affine Cartan matrix.
func Pair(a, b Vec) int64 {
	var s int64
	for i := 0; i < N; i++ {
		inner := 2 * b[i]
		for _, j := range Adj[i] {
			inner -= b[j]
		}
		Assert(abs(a[i]) <= limit && abs(inner) <= limit, "coefficient bound exceeded")
		s += a[i] * inner
	}
	return s
}

// MV is the matrix-vector product (rows dotted with v).
func MV(a Mat, v Vec) Vec {
	var out Vec
	for i := range out {
		out[i] = Dot(Vec(a[i]), v)
	}
	return out
}

// MM is the matrix product.
func MM(a, b Mat) Mat {
	var out Mat
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			var s int64
			for k := 0; k < N; k++ {
				Assert(abs(a[i][k]) <= limit && abs(b[k][j]) <= limit, "coefficient bound exceeded")
				s += a[i][k] * b[k][j]
			}
			out[i][j] = s
		}
	}
	return out
}

// Col returns column j of a.
func Col(a Mat, j int) Vec {
	var out Vec
	for i := range out {
		out[i] = a[i][j]
	}
	return out
}

// Positive asserts that a is a nonzero vector with coefficients of one sign
// and reports whether they are nonnegative.
func Positive(a Vec) bool {
	anyNonzero, allNonneg, allNonpos := false, true, true
	for _, c := range a {
		if c != 0 {
			anyNonzero = true
		}
		if c < 0 {
			allNonneg = false
		}
		if c > 0 {
			allNonpos = false
		}
	}
	Assert(anyNonzero, "zero vector %v", a)
	Assert(allNonneg || allNonpos, "mixed signs %v", a)
	return allNonneg
}

// Right is right multiplication by the simple reflection s.
func Right(a Mat, s int) Mat {
	out := a
	for i := 0; i < N; i++ {
		out[i][s] = -a[i][s]
		for _, t := range Adj[s] {
			out[i][t] = a[i][t] + a[i][s]
		}
	}
	for i := range out {
		for j := range out[i] {
			Assert(abs(out[i][j]) <= limit, "coefficient bound exceeded")
		}
	}
	return out
}

// Left is left multiplication by the simple reflection s (acts on row s).
func Left(a Mat, s int) Mat {
	out := a
	for j := 0; j < N; j++ {
		v := -a[s][j]
		for _, t := range Adj[s] {
			v += a[t][j]
		}
		out[s][j] = v
	}
	return out
}

// WordMatrix multiplies the generators of word on the right; with reduced
// set it asserts that every step increases the length.
func WordMatrix(word []int, reduced bool) Mat {
	a := E
	for _, s := range word {
		if reduced {
			Assert(Positive(Col(a, s)), "nonreduced %v", word)
		}
		a = Right(a, s)
	}
	return a
}

// ReducedWord returns the reduced word of a (lowest-index right descent first).
func ReducedWord(a Mat) []int {
	rev := []int{}
	steps := 0
	for a != E {
		steps++
		Assert(steps <= 100000, "reduced_word does not terminate")
		choice := -1
		for s := 0; s < N; s++ {
			if !Positive(Col(a, s)) {
				choice = s
				break
			}
		}
		Assert(choice >= 0, "no descent")
		rev = append(rev, choice)
		a = Right(a, choice)
	}
	out := make([]int, len(rev))
	for i, s := range rev {
		out[len(rev)-1-i] = s
	}
	return out
}

// Reflection is the reflection in the root a.
func Reflection(a Vec) Mat {
	var p Vec
	for j := 0; j < N; j++ {
		p[j] = Pair(Vec(E[j]), a)
	}
	var out Mat
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			out[i][j] = E[i][j] - a[i]*p[j]
		}
	}
	return out
}

// Translate is t_gamma^k with t_gamma(alpha)=alpha-<alpha,gamma>*delta.
func Translate(gamma Vec, k int64) Mat {
	var p Vec
	for j := 0; j < N; j++ {
		p[j] = Pair(Vec(E[j]), gamma)
	}
	var out Mat
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			out[i][j] = E[i][j] - k*Delta[i]*p[j]
		}
	}
	return out
}

// Shifted is a + k * delta * d (outer product).
func Shifted(a Mat, d Vec, k int64) Mat {
	var out Mat
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			out[i][j] = a[i][j] + k*Delta[i]*d[j]
		}
	}
	return out
}

// FiniteRoots is the set of the 240 finite E8 roots (node 7 removed) and
// FiniteRootList the same roots in a fixed sorted order.
var FiniteRoots, FiniteRootList = finiteRoots()

func vecLess(a, b Vec) bool {
	for i := range a {
		if a[i] != b[i] {
			return a[i] < b[i]
		}
	}
	return false
}

func finiteRoots() (map[Vec]bool, []Vec) {
	roots := map[Vec]bool{}
	todo := []Vec{}
	for i := 0; i < N; i++ {
		if i != 7 {
			v := Vec(E[i])
			roots[v] = true
			todo = append(todo, v)
		}
	}
	for idx := 0; idx < len(todo); idx++ {
		a := todo[idx]
		for s := 0; s < N; s++ {
			if s == 7 {
				continue
			}
			b := a
			b[s] = -a[s]
			for _, t := range Adj[s] {
				b[s] += a[t]
			}
			if !roots[b] {
				roots[b] = true
				todo = append(todo, b)
			}
		}
	}
	Assert(len(roots) == 240, "finite root count %d", len(roots))
	list := make([]Vec, 0, len(roots))
	for a := range roots {
		Assert(a[7] == 0 && Pair(a, a) == 2, "bad finite root %v", a)
		list = append(list, a)
	}
	sort.Slice(list, func(i, j int) bool { return vecLess(list[i], list[j]) })
	return roots, list
}

// InversionClass is one (A,B) class of finite roots with its multiplicity.
type InversionClass struct {
	A     int64 `json:"A"`
	B     int64 `json:"B"`
	Count int64 `json:"count"`
}

// LengthFormula proves ell(a+k delta*d)=intercept+slope*k for all k>=0.
type LengthFormula struct {
	Intercept            int64            `json:"intercept"`
	Slope                int64            `json:"slope"`
	InversionClasses     []InversionClass `json:"inversion_classes"`
	NoBreakpointsForKGe0 bool             `json:"no_breakpoints_for_k_ge_0"`
}

func b2i(b bool) int64 {
	if b {
		return 1
	}
	return 0
}

// ComputeLengthFormula ports length_formula.
func ComputeLengthFormula(a Mat, d Vec) LengthFormula {
	type key struct{ A, B int64 }
	classes := map[key]int64{}
	var intercept, slope int64
	for _, r := range FiniteRootList {
		image := MV(a, r)
		c := image[7]
		finiteImage := Add(image, Delta, -c)
		Assert(FiniteRoots[finiteImage], "image not a finite root")
		nmin := b2i(!Positive(r))
		A := -c - nmin + b2i(!Positive(finiteImage))
		B := -Dot(d, r)
		Assert(B <= 0 || A >= 0, "breakpoint %d %d %v", A, B, r)
		Assert(B >= 0 || A <= 0, "breakpoint %d %d %v", A, B, r)
		if A > 0 {
			intercept += A
		}
		if B > 0 {
			slope += B
		}
		classes[key{A, B}]++
	}
	keys := make([]key, 0, len(classes))
	for k := range classes {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		if keys[i].A != keys[j].A {
			return keys[i].A < keys[j].A
		}
		return keys[i].B < keys[j].B
	})
	out := LengthFormula{Intercept: intercept, Slope: slope, NoBreakpointsForKGe0: true,
		InversionClasses: make([]InversionClass, 0, len(keys))}
	for _, k := range keys {
		out.InversionClasses = append(out.InversionClasses, InversionClass{k.A, k.B, classes[k]})
	}
	return out
}

// TerminalEdge is one terminal edge check.
type TerminalEdge struct {
	S        int   `json:"s"`
	T        int   `json:"t"`
	Slope    int64 `json:"slope"`
	BaseRoot Vec   `json:"base_root"`
}

// TerminalCert is the symbolic terminal certificate.
type TerminalCert struct {
	Descents           []int          `json:"descents"`
	InvolutionsForAllK bool           `json:"involutions_for_all_k"`
	TerminalForAllK    bool           `json:"terminal_for_all_k"`
	TerminalEdgeChecks []TerminalEdge `json:"terminal_edge_checks"`
}

// ComputeTerminalCertificate ports terminal_certificate.
func ComputeTerminalCertificate(a Mat, d Vec) TerminalCert {
	Assert(MM(a, a) == E && MV(a, Delta) == Delta, "a is not an involution fixing delta")
	Assert(Dot(d, Delta) == 0, "d.delta != 0")
	for j := 0; j < N; j++ {
		var s int64
		for i := 0; i < N; i++ {
			s += d[i] * a[i][j]
		}
		Assert(s == -d[j], "d a != -d at column %d", j)
	}
	cert := TerminalCert{Descents: []int{}, InvolutionsForAllK: true, TerminalForAllK: true,
		TerminalEdgeChecks: []TerminalEdge{}}
	for s := 0; s < N; s++ {
		if Positive(Col(a, s)) {
			Assert(d[s] >= 0, "d[%d] < 0 at ascent", s)
		} else {
			cert.Descents = append(cert.Descents, s)
			Assert(d[s] <= 0, "d[%d] > 0 at descent", s)
		}
	}
	for _, s := range cert.Descents {
		for _, t := range Adj[s] {
			base := Add(Col(a, s), Col(a, t), 1)
			Assert(Positive(base), "base root not positive")
			Assert(d[s]+d[t] >= 0, "negative terminal edge slope")
			cert.TerminalEdgeChecks = append(cert.TerminalEdgeChecks,
				TerminalEdge{S: s, T: t, Slope: d[s] + d[t], BaseRoot: base})
		}
	}
	return cert
}

// Row is one element of research/en_independent/e9-fc.json.
type Row struct {
	Word   []int `json:"word"`
	Length int   `json:"length"`
	Rmask  int   `json:"Rmask"`
	Lmask  int   `json:"Lmask"`
}

// Catalogue is research/en_independent/e9-fc.json.
type Catalogue struct {
	Complete           bool    `json:"complete"`
	FcCount            int     `json:"fc_count"`
	AscentsTested      int64   `json:"ascents_tested"`
	LengthDistribution []int64 `json:"length_distribution"`
	Elements           []Row   `json:"elements"`
}

// CataloguePath is relative to the working directory.
const CataloguePath = "research/en_independent/e9-fc.json"

// SeedPath is the payload seed with the length-27 base words.
const SeedPath = "research/en_families/e9_full_support_eligible.json"

// LoadCatalogue reads the E9 FC catalogue.
func LoadCatalogue() Catalogue {
	raw, err := os.ReadFile(CataloguePath)
	if err != nil {
		panic(err)
	}
	var c Catalogue
	if err := json.Unmarshal(raw, &c); err != nil {
		panic(err)
	}
	return c
}

// BaseWord27 returns the first word of length 27 in the seed file.
func BaseWord27() []int {
	raw, err := os.ReadFile(SeedPath)
	if err != nil {
		panic(err)
	}
	var rows []struct {
		Word   []int `json:"word"`
		Length int   `json:"length"`
	}
	if err := json.Unmarshal(raw, &rows); err != nil {
		panic(err)
	}
	for _, r := range rows {
		if r.Length == 27 {
			return r.Word
		}
	}
	panic(AssertionError{"no length-27 word in seed"})
}

// DescentBraidWitness is a reduced expression for s*c_0 with an explicit
// noncommuting braid. BraidStart is a zero-based index into Word.
type DescentBraidWitness struct {
	Generator  int   `json:"generator"`
	Word       []int `json:"word"`
	BraidStart int   `json:"braid_start_zero_based"`
}

// Family is the audit result for one affine family.
type Family struct {
	Name                          string        `json:"name"`
	BaseWord                      []int         `json:"base_word"`
	BaseRowMatrix                 Mat           `json:"base_row_matrix"`
	ColumnSlopes                  Vec           `json:"column_slopes"`
	LengthFormula                 LengthFormula `json:"length_formula"`
	RightTranslationReducedWord   []int         `json:"right_translation_reduced_word"`
	RightTranslationLength        int           `json:"right_translation_length"`
	RightTranslationLengthFormula LengthFormula `json:"right_translation_length_formula"`
	RightTranslationPowersForAllK bool          `json:"right_translation_powers_for_all_k"`
	BasePrefixForAllK             bool          `json:"base_prefix_for_all_k"`
	DistinctAndFullSupportForAllK bool          `json:"distinct_and_full_support_for_all_k"`
	TerminalCert
	MaximumIndependentSize    int                   `json:"maximum_independent_size"`
	LeftDescentBraidWitnesses []DescentBraidWitness `json:"left_descent_braid_witnesses"`
	DescentPrefixesForAllK    bool                  `json:"descent_prefixes_for_all_k"`
	NoFCCoversForAllK         bool                  `json:"no_FC_covers_for_all_k"`
}

func intsEqual(a, b []int) bool {
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

// VerifyDescentBraids checks a complete set of reduced braid witnesses for
// s*a, s in I. It uses exact matrix products and root signs, not FC membership.
func VerifyDescentBraids(a Mat, baseLength int, I []int, witnesses []DescentBraidWitness) {
	Assert(len(witnesses) == len(I), "descent braid witness count")
	for i, witness := range witnesses {
		Assert(witness.Generator == I[i], "descent braid generator %d != %d", witness.Generator, I[i])
		Assert(len(witness.Word) == baseLength-1, "descent braid word length")
		for _, s := range witness.Word {
			Assert(s >= 0 && s < N, "invalid braid word generator %d", s)
		}
		Assert(WordMatrix(witness.Word, true) == Left(a, witness.Generator), "descent braid product for %d", witness.Generator)
		p := witness.BraidStart
		Assert(p >= 0 && p+2 < len(witness.Word), "invalid braid position")
		s, t, u := witness.Word[p], witness.Word[p+1], witness.Word[p+2]
		adjacent := false
		for _, v := range Adj[s] {
			adjacent = adjacent || v == t
		}
		Assert(s == u && adjacent, "not a noncommuting braid at %d", p)
	}
}

// MaximumIndependentSize enumerates the 2^9 vertex subsets of the fixed E9
// diagram. This small graph check has no group-element or FC-catalogue input.
func MaximumIndependentSize() int {
	maximum := 0
	for mask := 0; mask < 1<<N; mask++ {
		independent := true
		for _, e := range Edges {
			if mask&(1<<e[0]) != 0 && mask&(1<<e[1]) != 0 {
				independent = false
				break
			}
		}
		if independent {
			size := 0
			for s := 0; s < N; s++ {
				if mask&(1<<s) != 0 {
					size++
				}
			}
			if size > maximum {
				maximum = size
			}
		}
	}
	return maximum
}

// AuditFamily certifies terminality, lengths and absence of FC covers for
// every k>=0 using reduced translation prefixes and descent braid witnesses.
// expected is (length intercept, slope). No FC catalogue is read.
func AuditFamily(baseWord []int, d Vec, expected [2]int64, I []int, witnesses []DescentBraidWitness, name string) Family {
	a := WordMatrix(baseWord, true)
	symbolic := ComputeTerminalCertificate(a, d)
	Assert(intsEqual(symbolic.Descents, I), "descents %v != %v", symbolic.Descents, I)
	formula := ComputeLengthFormula(a, d)
	Assert(formula.Intercept == expected[0] && formula.Slope == expected[1], "length formula %d,%d", formula.Intercept, formula.Slope)
	Assert(int64(len(baseWord)) == expected[0], "base word length")
	support := map[int]bool{}
	for _, s := range baseWord {
		support[s] = true
	}
	Assert(len(support) == N, "base word does not have full support")
	R := Shifted(E, d, 1)
	Assert(MM(a, R) == Shifted(a, d, 1), "aR != shifted a")
	rword := ReducedWord(R)
	Assert(int64(len(rword)) == expected[1] && WordMatrix(rword, true) == R, "right translation word")
	rformula := ComputeLengthFormula(E, d)
	Assert(rformula.Intercept == 0 && rformula.Slope == expected[1], "right translation length formula")
	// d.delta=0, already checked by ComputeTerminalCertificate, gives
	// (R-1)^2=0 and R^k=1+k*delta*d on the whole affine root space.
	Assert(MM(R, R) == Shifted(E, d, 2), "right translation square")
	Assert(expected[1] > 0, "family is not pairwise distinct")
	// Since a.delta=delta, aR^k=a+k*delta*d. The length formulas above
	// give ell(aR^k)=ell(a)+ell(R^k). For s in the constant descent set,
	// ell(saR^k)=ell(sa)+ell(R^k), so the checked braid in sa persists.
	VerifyDescentBraids(a, len(baseWord), I, witnesses)
	maximum := MaximumIndependentSize()
	Assert(len(I) == maximum && len(I)%2 == 1, "descents not maximum of odd size")
	mask := 0
	for _, s := range I {
		mask |= 1 << s
	}
	for _, e := range Edges {
		Assert(mask&(1<<e[0]) == 0 || mask&(1<<e[1]) == 0, "descents not independent")
	}
	Assert(expected[0]%2 == int64(len(I)%2) && expected[1]%2 == 0, "maximum-descent parity")
	Assert(expected[0]-1 > int64(len(I)), "commuting product could be a base cover")
	// The maximum-descent lemma excludes an FC cover retaining every
	// left/right descent. The lifting property makes any other cover saR^k
	// or its inverse, both excluded by the persistent braid witnesses.
	return Family{Name: name, BaseWord: baseWord, BaseRowMatrix: a, ColumnSlopes: d,
		LengthFormula: formula, RightTranslationReducedWord: rword, RightTranslationLength: len(rword),
		RightTranslationLengthFormula: rformula, RightTranslationPowersForAllK: true,
		BasePrefixForAllK: true, DistinctAndFullSupportForAllK: true, TerminalCert: symbolic,
		MaximumIndependentSize: maximum, LeftDescentBraidWitnesses: witnesses,
		DescentPrefixesForAllK: true, NoFCCoversForAllK: true}
}

// WriteJSON writes v as indented JSON (2 spaces, no HTML escaping, trailing
// newline) to path, or to stdout when path is empty.
func WriteJSON(path string, v any) {
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(v); err != nil {
		panic(err)
	}
	if path == "" {
		if _, err := os.Stdout.Write(buf.Bytes()); err != nil {
			panic(err)
		}
		return
	}
	// Written in one checked call so that a write error cannot leave a
	// truncated certificate behind a successful exit.
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		panic(err)
	}
}
