// Command terminals-flat enumerates, exactly, the terminal elements of the
// finite Weyl groups E6, E7 and E8 by one-sided parabolic pruning: minimal
// right cosets of the parabolic subgroup W_J (J = all generators but one)
// times the one-sided terminal elements of W_J. No full group enumeration
// beyond W_J is performed. It then separates the commuting terminals (fully
// commutative elements) from the non-commuting ones and lists, for each of
// the latter, the fully commutative "bottoms" below it.
//
// It ports research/en_e8/parabolic_terminals.cpp. Usage:
//
//	terminals-flat -rank N [-cosets K] > eN-validation.json   (N = 6, 7, 8)
//
// The optional -cosets K limits the number of cosets processed (the result is
// then marked "complete": false). The certificate is written to stdout; the
// timing field of the original is dropped.
//
// Imports: the standard library and go/internal/parabolic (shared with
// terminals-recursive, as in the original). It imports none of the other
// engines of the supplement.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"

	"github.com/tygern/kl/supplement/internal/parabolic"
)

type certificate struct {
	Type                      string          `json:"type"`
	Complete                  bool            `json:"complete"`
	Method                    string          `json:"method"`
	ParabolicOmittedGenerator int             `json:"parabolic_omitted_generator"`
	ParabolicOrder            uint64          `json:"parabolic_order"`
	ParabolicRightTerminals   int             `json:"parabolic_right_terminal_count"`
	CosetCount                int             `json:"coset_count"`
	CosetsProcessed           int             `json:"cosets_processed"`
	CandidatesTested          uint64          `json:"candidates_tested"`
	RightTerminalCount        uint64          `json:"right_terminal_count"`
	FCCount                   int             `json:"fc_count"`
	FCStarChecks              uint64          `json:"fc_star_checks"`
	CommutingTerminals        int             `json:"commuting_terminals"`
	NoncommutingTerminalCount int             `json:"noncommuting_terminal_count"`
	Bad                       []parabolic.Bad `json:"bad"`
}

func main() {
	defer parabolic.ExitOnFailure()
	rank := flag.Int("rank", 8, "rank n of E_n (6, 7 or 8)")
	cosets := flag.Int("cosets", -1, "limit the number of cosets processed (default: all)")
	flag.Parse()
	parabolic.Assert(flag.NArg() == 0, "unexpected positional arguments")
	g := parabolic.New(*rank)
	n := g.N
	par, parSize := g.ParabolicOneSided()
	reps := g.Cosets()
	fc := g.FullyCommutative()
	limit := len(reps)
	if *cosets >= 0 && *cosets < limit {
		limit = *cosets
	}
	all := uint32(1)<<uint(n) - 1
	var candidates, rightTerminals uint64
	var terminals []parabolic.Element
	for ci := 0; ci < limit; ci++ {
		a := reps[ci]
		perm := g.RootPermutation(a)
		for _, v := range par {
			candidates++
			w := g.Compose(perm, v.W)
			if !g.WeakRight(w, all) {
				continue
			}
			rightTerminals++
			wi := v.Inv
			for i := len(a.Word) - 1; i >= 0; i-- {
				wi = g.Right(wi, a.Word[i])
			}
			if !g.WeakRight(wi, all) {
				continue
			}
			word := append(append([]int(nil), a.Word...), v.Word...)
			parabolic.Assert(len(g.ReducedWord(w)) == len(word), "word is not reduced")
			parabolic.Assert(g.Inverse(word) == wi, "inverse mismatch")
			terminals = append(terminals, parabolic.Element{W: w, Inv: wi, Word: word})
		}
		if ci%20 == 19 || ci+1 == limit {
			fmt.Fprintf(os.Stderr, "Cosets %d/%d, candidates %d, right terminals %d, terminals %d\n",
				ci+1, limit, candidates, rightTerminals, len(terminals))
		}
	}
	parabolic.SortByLengthThenWord(terminals)
	fcIDs := parabolic.FCIndex(fc)
	commuting := 0
	var bads []parabolic.Element
	for _, b := range terminals {
		if _, ok := fcIDs[b.W]; ok {
			g.SupportCommuting(b.Word)
			commuting++
		} else {
			bads = append(bads, b)
		}
	}
	cert := certificate{
		Type:                      fmt.Sprintf("E%d", n),
		Complete:                  limit == len(reps),
		Method:                    "minimal right parabolic cosets times one-sided terminal parabolic elements",
		ParabolicOmittedGenerator: g.Omitted,
		ParabolicOrder:            parSize,
		ParabolicRightTerminals:   len(par),
		CosetCount:                len(reps),
		CosetsProcessed:           limit,
		CandidatesTested:          candidates,
		RightTerminalCount:        rightTerminals,
		FCCount:                   len(fc),
		FCStarChecks:              g.FCStarChecks,
		CommutingTerminals:        commuting,
		NoncommutingTerminalCount: len(bads),
		Bad:                       g.BadRecords(bads, fc),
	}
	out, err := json.Marshal(cert)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	fmt.Println(string(out))
}
