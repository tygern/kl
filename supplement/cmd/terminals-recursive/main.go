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
// It ports research/en_e8/recursive_terminals.cpp. Usage:
//
//	terminals-recursive -rank N > eN-recursive.json   (N = 6, 7, 8)
//
// The certificate is written to stdout; the timing and peak-memory fields of
// the original are dropped.
//
// Imports: the standard library and go/internal/parabolic (shared with
// terminals-flat, as in the original, which #includes the flat program as a
// library). It imports none of the other engines of the supplement.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"

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

func main() {
	defer parabolic.ExitOnFailure()
	rank := flag.Int("rank", 8, "rank n of E_n (6, 7 or 8)")
	flag.Parse()
	parabolic.Assert(flag.NArg() == 0, "unexpected positional arguments")
	g := parabolic.New(*rank)
	n := g.N
	additions := []int{n - 1, 2, 3, 1, 0, 4}
	for s := 5; s <= n-2; s++ {
		additions = append(additions, s)
	}
	indices := []int64{2, 3, 4, 8, 10, 27, 56, 240}
	parabolic.Assert(len(additions) == n, "addition chain has wrong length")
	rightTerminals := []parabolic.Element{{W: g.Identity, Inv: g.Identity, Word: nil}}
	stages := []stage{}
	var K uint32
	var order, totalCandidates int64 = 1, 0
	for level := range additions {
		H := K
		K |= 1 << uint(additions[level])
		previous := len(rightTerminals)
		next := extendRightTerminals(g, rightTerminals, K, H, int(indices[level]))
		order = parabolic.MulChecked(order, indices[level])
		cands := parabolic.MulChecked(indices[level], int64(previous))
		totalCandidates += cands
		parabolic.Assert(totalCandidates >= 0, "candidate count overflow")
		stages = append(stages, stage{K, H, additions[level], order, indices[level], previous, cands, len(next)})
		fmt.Fprintf(os.Stderr, "Rank %d, order %d, cosets %d, candidates %d, right terminals %d\n",
			level+1, order, indices[level], cands, len(next))
		rightTerminals = next
	}
	wantOrder := int64(696729600)
	if n == 6 {
		wantOrder = 51840
	} else if n == 7 {
		wantOrder = 2903040
	}
	parabolic.Assert(order == wantOrder, "group order %d, expected %d", order, wantOrder)
	fc := g.FullyCommutative()
	fcIDs := parabolic.FCIndex(fc)
	commuting := 0
	var bads []parabolic.Element
	for _, b := range rightTerminals {
		if !g.WeakRight(b.Inv, K) {
			continue
		}
		if _, ok := fcIDs[b.W]; ok {
			g.SupportCommuting(b.Word)
			commuting++
		} else {
			bads = append(bads, b)
		}
	}
	parabolic.SortByLengthThenWord(bads)
	cert := certificate{
		Type:               fmt.Sprintf("E%d", n),
		Complete:           true,
		Method:             "recursive one-sided parabolic pruning",
		GroupOrder:         order,
		TotalCandidates:    totalCandidates,
		RightTerminalCount: len(rightTerminals),
		FCCount:            len(fc),
		FCStarChecks:       g.FCStarChecks,
		CommutingTerminals: commuting,
		NoncommutingCount:  len(bads),
		Stages:             stages,
		Bad:                g.BadRecords(bads, fc),
	}
	out, err := json.Marshal(cert)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	fmt.Println(string(out))
}
