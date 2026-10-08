// Command e10-all-k is the independent all-k audit of the full-support E10
// reflection family (beta_k = T^k beta_0 with T-1 nilpotent of order 3).
//
// It ports research/en_affine_referee/verify_indefinite_e10.py, which in turn
// uses the library part of research/en_affine_referee/verify_indefinite.py;
// the shared library is go/internal/indefinite. It reads the payload seed
// research/en_families/cartan_E10_m2_max1.json relative to the working
// directory and writes research/en_affine_referee/e10-certificate.json.
// It imports no other engine of this module.
package main

import (
	"encoding/json"
	"errors"
	"fmt"
	"os"

	q "github.com/tygern/kl/supplement/internal/indefinite"
)

var (
	beta  = q.Vec{3, 7, 10, 9, 7, 6, 4, 3, 1, 6}
	delta = q.Vec{2, 4, 6, 5, 4, 3, 2, 1, 0, 3}
	gamma = q.Vec{1, 2, 3, 2, 2, 1, 1, 0, 0, 1}
	setI  = []int{1, 3, 5, 7, 9}
)

const (
	seedPath = "research/en_families/cartan_E10_m2_max1.json"
	outPath  = "research/en_affine_referee/e10-certificate.json"
)

type edgeCheck struct {
	S                int     `json:"s"`
	T                int     `json:"t"`
	PairingSumCoeffs []int64 `json:"pairing_sum_coefficients"`
}

type seedRow struct {
	Beta        []int64 `json:"beta"`
	ReducedWord []int   `json:"reduced_word"`
}

type seedFile struct {
	RealTerminalRoots []seedRow `json:"real_terminal_roots"`
}

type certificate struct {
	Status                 string               `json:"status"`
	Rank                   int                  `json:"rank"`
	Beta0                  q.Vec                `json:"beta0"`
	EmbeddedAffineDelta    q.Vec                `json:"embedded_affine_delta"`
	Gamma                  q.Vec                `json:"gamma"`
	RealRootWitnesses      map[string]q.Witness `json:"real_root_witnesses"`
	BaseLength             int                  `json:"base_length"`
	TMinusIdentityCubeZero bool                 `json:"T_minus_identity_cube_zero"`
	NBeta0                 q.Vec                `json:"N_beta0"`
	NSquaredBeta0          q.Vec                `json:"N_squared_beta0"`
	NUZero                 bool                 `json:"N_u_zero"`
	BetaFormula            string               `json:"beta_formula"`
	PairingCoefficients    []q.Vec              `json:"pairing_coefficients"`
	FixedLeftRightDescents []int                `json:"fixed_left_right_descents"`
	TerminalEdgeChecks     []edgeCheck          `json:"terminal_edge_checks"`
	MaximumIndependentSize int                  `json:"maximum_independent_size"`
	FullSupportForAllK     bool                 `json:"full_support_for_all_k"`
	FullSupportWitness     string               `json:"full_support_witness"`
	DistinctnessWitness    string               `json:"distinctness_witness"`
	LengthParity           string               `json:"length_parity"`
	UniqueEligibleFCBottom []int                `json:"unique_eligible_FC_bottom"`
	EligibleMu             int                  `json:"eligible_mu"`
	AllFCLowerMuBound      []int                `json:"all_FC_lower_mu_bound"`
	ExactLengthFormula     string               `json:"exact_length_formula"`
}

func inI(s int) bool {
	for _, x := range setI {
		if x == s {
			return true
		}
	}
	return false
}

func equal(a q.Vec, b ...int64) bool { return q.VecEq(a, q.Vec(b)) }

func run() error {
	c := q.Configure(10)
	N := c.N
	witnesses := map[string]q.Witness{
		"beta0":            c.RealRootWitness(beta),
		"gamma":            c.RealRootWitness(gamma),
		"gamma_plus_delta": c.RealRootWitness(q.Add(gamma, delta, 1)),
	}
	raw, err := os.ReadFile(seedPath)
	if err != nil {
		return err
	}
	var source seedFile
	if err := json.Unmarshal(raw, &source); err != nil {
		return err
	}
	var row *seedRow
	for i := range source.RealTerminalRoots {
		if q.VecEq(source.RealTerminalRoots[i].Beta, beta) {
			row = &source.RealTerminalRoots[i]
			break
		}
	}
	if row == nil {
		return errors.New("seed has no real terminal root equal to beta0")
	}
	q.Assert(len(row.ReducedWord) == 101, "reduced word length is %d, not 101", len(row.ReducedWord))
	q.Assert(q.MatEq(c.ReducedMatrix(row.ReducedWord), c.Reflection(beta)), "reduced word does not give the reflection")
	letters := map[int]bool{}
	for _, s := range row.ReducedWord {
		letters[s] = true
	}
	allLetters := len(letters) == N
	for s := 0; s < N; s++ {
		allLetters = allLetters && letters[s]
	}
	q.Assert(allLetters, "reduced word does not use every generator")

	T := c.MM(c.Reflection(gamma), c.Reflection(q.Add(gamma, delta, 1)))
	M := make(q.Mat, N)
	for i := 0; i < N; i++ {
		M[i] = make([]int64, N)
		for j := 0; j < N; j++ {
			M[i][j] = T[i][j] - c.E[i][j]
		}
	}
	M2 := c.MM(M, M)
	q.Assert(q.MatEq(c.MM(M2, M), c.Zero), "(T-1)^3 is not zero")
	v := c.MV(M, beta)
	u := c.MV(M2, beta)
	q.Assert(equal(v, 7, 14, 21, 18, 14, 11, 7, 4, 0, 11), "N beta0 mismatch")
	q.Assert(q.VecEq(u, q.Add(make(q.Vec, N), delta, 2)), "N^2 beta0 != 2 delta")
	q.Assert(equal(c.MV(M, u), make(q.Vec, N)...), "N u != 0")
	q.Assert(q.Min(beta) > 0 && q.Min(v) >= 0 && q.Min(u) >= 0, "positivity failed")
	m0, m1, m2 := c.PairingVector(beta), c.PairingVector(v), c.PairingVector(u)
	q.Assert(equal(m0, -1, 1, -2, 1, -1, 1, -1, 1, -1, 2), "m0 mismatch")
	q.Assert(equal(m1, 0, 0, -1, 1, -1, 1, -1, 1, -4, 1), "m1 mismatch")
	q.Assert(equal(m2, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0), "m2 mismatch")
	for s := 0; s < N; s++ {
		co := []int64{m0[s], m1[s], m2[s]}
		if inI(s) {
			q.Assert(co[0] >= 1 && co[1] >= 0 && co[2] >= 0, "coefficient signs at s=%d in I", s)
		} else {
			q.Assert(co[0] <= 0 && co[1] <= 0 && co[2] <= 0, "coefficient signs at s=%d not in I", s)
		}
	}
	edgeChecks := []edgeCheck{}
	for _, s := range setI {
		for _, t := range c.Adj[s] { // sorted ascending
			co := []int64{m0[s] + m0[t], m1[s] + m1[t], m2[s] + m2[t]}
			q.Assert(co[0] <= 0 && co[1] <= 0 && co[2] <= 0, "edge pairing sum positive at (%d,%d)", s, t)
			edgeChecks = append(edgeChecks, edgeCheck{S: s, T: t, PairingSumCoeffs: co})
		}
	}
	maxIndep := 0
	for m := 0; m < 1<<N; m++ {
		ok := true
		for _, e := range c.Edges {
			if m&(1<<e[0]) != 0 && m&(1<<e[1]) != 0 {
				ok = false
				break
			}
		}
		if ok {
			pc := 0
			for x := m; x != 0; x &= x - 1 {
				pc++
			}
			if pc > maxIndep {
				maxIndep = pc
			}
		}
	}
	q.Assert(maxIndep == len(setI) && len(setI) == 5, "maximum independent size is %d", maxIndep)
	for _, e := range c.Edges {
		q.Assert(!(inI(e[0]) && inI(e[1])), "I is not independent")
	}
	q.Assert(m0[1] == 1 && m1[1] == 0 && m2[1] == 0, "pairing with alpha1 mismatch")
	q.Assert(beta[0] == 3 && v[0] == 7 && u[0] == 4, "first coordinates mismatch")

	out := certificate{
		Status: "All exact all-k assertions passed", Rank: N,
		Beta0: beta, EmbeddedAffineDelta: delta, Gamma: gamma,
		RealRootWitnesses: witnesses, BaseLength: 101,
		TMinusIdentityCubeZero: true, NBeta0: v, NSquaredBeta0: u, NUZero: true,
		BetaFormula:            "beta_k=beta0+k*v+binom(k,2)*u=T^k*beta0",
		PairingCoefficients:    []q.Vec{m0, m1, m2},
		FixedLeftRightDescents: setI, TerminalEdgeChecks: edgeChecks,
		MaximumIndependentSize: 5, FullSupportForAllK: true,
		FullSupportWitness:     "<alpha1,beta_k>=1, so column1 of r_beta_k differs from identity in every row",
		DistinctnessWitness:    "beta_k[0]=2*k*k+5*k+3",
		LengthParity:           "odd for every k, because each element is a real-root reflection",
		UniqueEligibleFCBottom: setI, EligibleMu: 0,
		AllFCLowerMuBound:  []int{0, 1},
		ExactLengthFormula: "not claimed or required",
	}
	data, err := json.MarshalIndent(out, "", "  ")
	if err != nil {
		return err
	}
	if err := os.WriteFile(outPath, append(data, '\n'), 0o644); err != nil {
		return err
	}
	summary := map[string]any{
		"status": out.Status, "rank": out.Rank, "full_support_for_all_k": out.FullSupportForAllK,
		"unique_eligible_FC_bottom": out.UniqueEligibleFCBottom, "eligible_mu": out.EligibleMu,
		"all_FC_lower_mu_bound": out.AllFCLowerMuBound,
	}
	sd, err := json.MarshalIndent(summary, "", "  ")
	if err != nil {
		return err
	}
	fmt.Println(string(sd))
	return nil
}

func main() {
	defer func() {
		if r := recover(); r != nil {
			fmt.Fprintf(os.Stderr, "e10-all-k: %v\n", r)
			os.Exit(1)
		}
	}()
	if err := run(); err != nil {
		fmt.Fprintf(os.Stderr, "e10-all-k: %v\n", err)
		os.Exit(1)
	}
}
