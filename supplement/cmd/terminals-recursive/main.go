// Command terminals-recursive enumerates, exactly, the terminal elements of
// the finite Weyl groups E6, E7 and E8 by recursive one-sided parabolic
// pruning along a chain of standard parabolic subgroups
// W_{K_1} < W_{K_2} < ... < W_{K_n} = W, adding one generator at a time.
// At each stage the right-terminal elements of the larger subgroup are the
// products (minimal coset representative) * (right-terminal element of the
// smaller subgroup) that are still right terminal; the full parabolic
// subgroup is never enumerated. The computation also checks the
// sign-preservation of minimal representatives used in the proof.
//
// It ports research/en_e8/recursive_terminals.cpp and, for the E8 chains,
// the pruning mode of the review program en.cpp (an independent exact-integer
// implementation that was run along four parabolic chains; it is in the
// repository history at commit d997526).
// Usage:
//
//	terminals-recursive -rank N > eN-recursive.json   (N = 6, 7, 8)
//	terminals-recursive -rank 8 -chain e7|d7|a7|alt > chain.json
//	terminals-recursive -rank 8 -chains-certificate research/en_e8/e8-recursive-chains.json
//
// The chains differ in the order in which generators are added (node
// numbering: path 0-1-...-(n-2), node n-1 attached to node 2):
//
//	e7   (default, the paper's chain)  7,2,3,1,0,4,5,6   last step E7 < E8
//	d7   6,5,4,3,2,1,7,0                                  last step D7 < E8
//	a7   0,1,2,3,4,5,6,7                                  last step A7 < E8
//	alt  0,1,7,2,3,4,5,6                                  last step E7 < E8
//
// With -chain the recursive certificate of that chain is written to stdout.
// The default chain e7 gives byte for byte the default output; the other
// chains add a "chain" key. With -chains-certificate all four chains are run
// and, per chain, the coset counts, right-terminal counts, candidate counts
// and the sorted canonical words of the two-sided terminal elements are
// written to the named file; the program asserts that the right-terminal
// sets and the two-sided terminal sets of all four chains are equal as group
// elements. The timing and peak-memory fields of the original are dropped.
//
// Coset counts are not tabulated: the index of W_H in W_K is computed from
// the orders of the two parabolic subgroups (classified from their diagrams),
// and the coset enumeration asserts that it finds exactly that many cosets.
//
// Imports: the standard library and internal/parabolic (shared with
// terminals-flat, as in the original, which #includes the flat program as a
// library). It imports none of the other engines of the supplement.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"sort"
	"strings"

	"github.com/tygern/kl/supplement/internal/parabolic"
)

type stage struct {
	K        uint32 `json:"subgroup_mask"`
	J        uint32 `json:"smaller_subgroup_mask"`
	Added    int    `json:"added_generator"`
	Order    int64  `json:"order"`
	Cosets   int64  `json:"coset_count"`
	Previous int    `json:"previous_right_terminal_count"`
	Cands    int64  `json:"candidates_tested"`
	RightCnt int    `json:"right_terminal_count"`
}

type certificate struct {
	Type               string          `json:"type"`
	Complete           bool            `json:"complete"`
	Method             string          `json:"method"`
	Chain              string          `json:"chain,omitempty"`
	GroupOrder         int64           `json:"group_order"`
	TotalCandidates    int64           `json:"total_candidates_tested"`
	RightTerminalCount int             `json:"right_terminal_count"`
	FCCount            int             `json:"fc_count"`
	FCStarChecks       uint64          `json:"fc_star_checks"`
	CommutingTerminals int             `json:"commuting_terminals"`
	NoncommutingCount  int             `json:"noncommuting_terminal_count"`
	Stages             []stage         `json:"stages"`
	Bad                []parabolic.Bad `json:"bad"`
}

// chainRecord is one chain's entry of the multi-chain certificate.
type chainRecord struct {
	Chain               string   `json:"chain"`
	AddedGenerators     []int    `json:"added_generators"`
	GroupOrder          int64    `json:"group_order"`
	CosetCounts         []int64  `json:"coset_counts"`
	RightTerminalCounts []int    `json:"right_terminal_counts"`
	Candidates          []int64  `json:"candidates_tested"`
	TotalCandidates     int64    `json:"total_candidates_tested"`
	RightTerminalCount  int      `json:"right_terminal_count"`
	TerminalCount       int      `json:"terminal_count"`
	CommutingTerminals  int      `json:"commuting_terminals"`
	NoncommutingCount   int      `json:"noncommuting_terminal_count"`
	TerminalWords       []string `json:"terminal_words"`
}

type chainsCertificate struct {
	Type                   string        `json:"type"`
	Complete               bool          `json:"complete"`
	Method                 string        `json:"method"`
	WordConvention         string        `json:"word_convention"`
	GroupOrder             int64         `json:"group_order"`
	RightTerminalCount     int           `json:"right_terminal_count"`
	TerminalCount          int           `json:"terminal_count"`
	CommutingTerminals     int           `json:"commuting_terminals"`
	NoncommutingCount      int           `json:"noncommuting_terminal_count"`
	RightTerminalSetsEqual bool          `json:"right_terminal_sets_equal"`
	TerminalSetsEqual      bool          `json:"terminal_sets_equal"`
	Chains                 []chainRecord `json:"chains"`
}

// chainRun is the outcome of the recursive pruning along one chain.
type chainRun struct {
	name            string
	additions       []int
	stages          []stage
	order           int64
	totalCandidates int64
	K               uint32
	rightTerminals  []parabolic.Element
}

// extendRightTerminals extends the right-terminal elements of W_H (prior) to
// those of W_K, where index is the number of minimal cosets of W_H in W_K.
func extendRightTerminals(g *parabolic.Group, prior []parabolic.Element, K, H uint32, index int) []parabolic.Element {
	reps := g.MinimalCosets(K, H, index)
	var result []parabolic.Element
	seen := map[parabolic.State]int{}
	for _, a := range reps {
		perm := g.RootPermutation(a)
		// Computational check of the sign-preservation used in the proof.
		for r := 0; r < g.NumRoots(); r++ {
			if !g.Positive(r) {
				continue
			}
			inSubsystem := true
			for s := 0; s < g.N; s++ {
				if H&(1<<uint(s)) == 0 && g.Roots[r][s] != 0 {
					inSubsystem = false
				}
			}
			if inSubsystem {
				parabolic.Assert(g.Positive(int(perm[r])), "sign preservation fails")
			}
		}
		for _, v := range prior {
			w := g.Compose(perm, v.W)
			// A minimal representative preserves all signs in the H subsystem.
			parabolic.Assert(g.WeakRight(w, H), "representative does not preserve H-terminality")
			if !g.WeakRight(w, K) {
				continue
			}
			_, dup := seen[w]
			parabolic.Assert(!dup, "duplicate product element")
			seen[w] = len(result)
			wi := v.Inv
			for i := len(a.Word) - 1; i >= 0; i-- {
				wi = g.Right(wi, a.Word[i])
			}
			word := append(append([]int(nil), a.Word...), v.Word...)
			parabolic.Assert(len(g.ReducedWord(w)) == len(word), "word is not reduced")
			parabolic.Assert(g.Inverse(word) == wi, "inverse mismatch")
			result = append(result, parabolic.Element{W: w, Inv: wi, Word: word})
		}
	}
	return result
}

// runChain prunes along the named chain, adding the generators of the chain
// one at a time, and asserts that the orders multiply to |W|.
func runChain(g *parabolic.Group, name string) chainRun {
	additions := g.AdditionChain(name)
	n := g.N
	parabolic.Assert(len(additions) == n, "addition chain has wrong length")
	var mask uint32
	for _, a := range additions {
		parabolic.Assert(a >= 0 && a < n && mask&(1<<uint(a)) == 0, "addition chain is not a permutation of the nodes")
		mask |= 1 << uint(a)
	}
	rightTerminals := []parabolic.Element{{W: g.Identity, Inv: g.Identity, Word: nil}}
	stages := []stage{}
	var K uint32
	var order, totalCandidates int64 = 1, 0
	for level := range additions {
		H := K
		K |= 1 << uint(additions[level])
		orderK, orderH := g.ParabolicOrder(K), g.ParabolicOrder(H)
		parabolic.Assert(orderK%orderH == 0, "parabolic order %d not divisible by %d", orderK, orderH)
		index := orderK / orderH
		previous := len(rightTerminals)
		next := extendRightTerminals(g, rightTerminals, K, H, int(index))
		order = parabolic.MulChecked(order, index)
		cands := parabolic.MulChecked(index, int64(previous))
		totalCandidates += cands
		parabolic.Assert(totalCandidates >= 0, "candidate count overflow")
		stages = append(stages, stage{K, H, additions[level], order, index, previous, cands, len(next)})
		fmt.Fprintf(os.Stderr, "Rank %d, order %d, cosets %d, candidates %d, right terminals %d\n",
			level+1, order, index, cands, len(next))
		rightTerminals = next
	}
	wantOrder := int64(696729600)
	if n == 6 {
		wantOrder = 51840
	} else if n == 7 {
		wantOrder = 2903040
	}
	parabolic.Assert(order == wantOrder, "group order %d, expected %d", order, wantOrder)
	return chainRun{name, additions, stages, order, totalCandidates, K, rightTerminals}
}

// terminalClassification splits the two-sided terminals of a run into the
// commuting products (asserting their support condition) and the
// noncommuting ones, the latter sorted by length then word.
type terminalClassification struct {
	twoSided   []parabolic.Element
	commuting  int
	noncommute []parabolic.Element
}

func classify(g *parabolic.Group, run chainRun, fcIDs map[parabolic.State]int) terminalClassification {
	var c terminalClassification
	for _, b := range run.rightTerminals {
		if !g.WeakRight(b.Inv, run.K) {
			continue
		}
		c.twoSided = append(c.twoSided, b)
		if _, ok := fcIDs[b.W]; ok {
			g.SupportCommuting(b.Word)
			c.commuting++
		} else {
			c.noncommute = append(c.noncommute, b)
		}
	}
	parabolic.SortByLengthThenWord(c.noncommute)
	return c
}

func writeJSON(v any, path string) {
	out, err := json.Marshal(v)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if path == "" {
		fmt.Println(string(out))
		return
	}
	if err := os.WriteFile(path, append(out, '\n'), 0o644); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

// sortedStates returns the sorted States of a list of elements.
func sortedStates(es []parabolic.Element) []parabolic.State {
	out := make([]parabolic.State, len(es))
	for i, e := range es {
		out[i] = e.W
	}
	sort.Slice(out, func(i, j int) bool { return out[i] < out[j] })
	for i := 1; i < len(out); i++ {
		parabolic.Assert(out[i] != out[i-1], "duplicate element in a terminal set")
	}
	return out
}

func sameStates(a, b []parabolic.State) bool {
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

// canonicalWord is the reduced word of w obtained by repeatedly stripping
// the lowest right descent, as a digit string (the node labels are 0..7).
func canonicalWord(g *parabolic.Group, w parabolic.State) string {
	var sb strings.Builder
	for _, s := range g.ReducedWord(w) {
		sb.WriteByte(byte('0' + s))
	}
	return sb.String()
}

func runAllChains(g *parabolic.Group, path string) {
	parabolic.Assert(g.N == 8, "-chains-certificate is defined for rank 8 only")
	fc := g.FullyCommutative()
	fcIDs := parabolic.FCIndex(fc)
	cert := chainsCertificate{
		Type:     "E8",
		Complete: true,
		Method:   "recursive one-sided parabolic pruning along four chains of parabolic subgroups",
		WordConvention: "terminal_words are the canonical reduced words of the two-sided terminal elements " +
			"(repeatedly strip the lowest right descent), sorted by length then lexicographically",
	}
	var refRight, refTerm []parabolic.State
	for ci, name := range parabolic.ChainNames {
		run := runChain(g, name)
		c := classify(g, run, fcIDs)
		right, term := sortedStates(run.rightTerminals), sortedStates(c.twoSided)
		if ci == 0 {
			refRight, refTerm = right, term
		}
		parabolic.Assert(sameStates(right, refRight), "right-terminal set of chain %s differs from chain %s", name, parabolic.ChainNames[0])
		parabolic.Assert(sameStates(term, refTerm), "terminal set of chain %s differs from chain %s", name, parabolic.ChainNames[0])
		words := make([]string, len(c.twoSided))
		for i, b := range c.twoSided {
			words[i] = canonicalWord(g, b.W)
		}
		sort.Slice(words, func(i, j int) bool {
			if len(words[i]) != len(words[j]) {
				return len(words[i]) < len(words[j])
			}
			return words[i] < words[j]
		})
		rec := chainRecord{
			Chain:               name,
			AddedGenerators:     run.additions,
			GroupOrder:          run.order,
			TotalCandidates:     run.totalCandidates,
			RightTerminalCount:  len(run.rightTerminals),
			TerminalCount:       len(c.twoSided),
			CommutingTerminals:  c.commuting,
			NoncommutingCount:   len(c.noncommute),
			TerminalWords:       words,
			CosetCounts:         []int64{},
			RightTerminalCounts: []int{},
			Candidates:          []int64{},
		}
		for _, st := range run.stages {
			rec.CosetCounts = append(rec.CosetCounts, st.Cosets)
			rec.RightTerminalCounts = append(rec.RightTerminalCounts, st.RightCnt)
			rec.Candidates = append(rec.Candidates, st.Cands)
		}
		// The values proved in the paper for E8 (Table 1 and the proposition
		// on the computational classification).
		parabolic.Assert(rec.RightTerminalCount == 2160 && rec.TerminalCount == 64 &&
			rec.CommutingTerminals == 58 && rec.NoncommutingCount == 6,
			"chain %s: counts %d/%d/%d/%d, expected 2160/64/58/6", name,
			rec.RightTerminalCount, rec.TerminalCount, rec.CommutingTerminals, rec.NoncommutingCount)
		cert.Chains = append(cert.Chains, rec)
		fmt.Fprintf(os.Stderr, "Chain %s: right terminals %d, terminals %d (commuting %d, noncommuting %d), candidates %d\n",
			name, rec.RightTerminalCount, rec.TerminalCount, rec.CommutingTerminals, rec.NoncommutingCount, rec.TotalCandidates)
	}
	first := cert.Chains[0]
	cert.GroupOrder = first.GroupOrder
	cert.RightTerminalCount = first.RightTerminalCount
	cert.TerminalCount = first.TerminalCount
	cert.CommutingTerminals = first.CommutingTerminals
	cert.NoncommutingCount = first.NoncommutingCount
	cert.RightTerminalSetsEqual = true
	cert.TerminalSetsEqual = true
	writeJSON(cert, path)
}

func main() {
	defer parabolic.ExitOnFailure()
	rank := flag.Int("rank", 8, "rank n of E_n (6, 7 or 8)")
	chain := flag.String("chain", "e7", "parabolic chain: e7 (default, the paper's), d7, a7 or alt (the last three for rank 8)")
	chainsPath := flag.String("chains-certificate", "", "run all four chains (rank 8) and write the multi-chain certificate to this file")
	flag.Parse()
	parabolic.Assert(flag.NArg() == 0, "unexpected positional arguments")
	g := parabolic.New(*rank)
	if *chainsPath != "" {
		parabolic.Assert(*chain == "e7", "-chain and -chains-certificate are mutually exclusive")
		runAllChains(g, *chainsPath)
		return
	}
	run := runChain(g, *chain)
	fc := g.FullyCommutative()
	c := classify(g, run, parabolic.FCIndex(fc))
	label := ""
	if *chain != "e7" {
		label = *chain
	}
	cert := certificate{
		Type:               fmt.Sprintf("E%d", g.N),
		Complete:           true,
		Method:             "recursive one-sided parabolic pruning",
		Chain:              label,
		GroupOrder:         run.order,
		TotalCandidates:    run.totalCandidates,
		RightTerminalCount: len(run.rightTerminals),
		FCCount:            len(fc),
		FCStarChecks:       g.FCStarChecks,
		CommutingTerminals: c.commuting,
		NoncommutingCount:  len(c.noncommute),
		Stages:             run.stages,
		Bad:                g.BadRecords(c.noncommute, fc),
	}
	writeJSON(cert, "")
}
