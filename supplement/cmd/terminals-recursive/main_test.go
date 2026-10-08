package main

import (
	"testing"

	"github.com/tygern/kl/supplement/internal/parabolic"
)

// TestChains prunes E8 along the four chains and checks the counts of the
// review and that the right-terminal sets coincide as group elements.
func TestChains(t *testing.T) {
	g := parabolic.New(8)
	want := map[string]struct {
		cosets     []int64
		candidates int64
	}{
		"e7":  {[]int64{2, 3, 4, 8, 10, 27, 56, 240}, 143660},
		"d7":  {[]int64{2, 3, 4, 5, 6, 7, 64, 2160}, 607290},
		"a7":  {[]int64{2, 3, 4, 5, 6, 7, 8, 17280}, 1210130},
		"alt": {[]int64{2, 3, 2, 10, 16, 27, 56, 240}, 143586},
	}
	var ref []parabolic.State
	for _, name := range parabolic.ChainNames {
		run := runChain(g, name)
		w := want[name]
		if run.totalCandidates != w.candidates || len(run.rightTerminals) != 2160 {
			t.Fatalf("chain %s: candidates %d, right terminals %d", name, run.totalCandidates, len(run.rightTerminals))
		}
		for i, st := range run.stages {
			if st.Cosets != w.cosets[i] {
				t.Fatalf("chain %s stage %d: %d cosets, want %d", name, i, st.Cosets, w.cosets[i])
			}
		}
		states := sortedStates(run.rightTerminals)
		if ref == nil {
			ref = states
		} else if !sameStates(ref, states) {
			t.Fatalf("chain %s: right-terminal set differs from chain e7", name)
		}
	}
}
