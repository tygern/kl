package main

import (
	"slices"
	"strings"
	"testing"
)

func TestCPythonHashes(t *testing.T) {
	if got := int64(intHash(-1)); got != -2 {
		t.Fatalf("hash(-1) = %d", got)
	}
	if got := int64(tupleHash([]uint64{1, 2, 3})); got != 529344067295497451 {
		t.Fatalf("hash((1,2,3)) = %d", got)
	}
	a := tupleHash([]uint64{intHash(-1), intHash(0)})
	b := tupleHash([]uint64{intHash(1), intHash(1)})
	if got := int64(tupleHash([]uint64{a, b})); got != -1287991002389685469 {
		t.Fatalf("nested tuple hash = %d", got)
	}
}

func TestSmallStars(t *testing.T) {
	// m = 3: P = 1 + 2q, |lower| = 4^3+2^3, |interval| = 3^3+1, w terminal.
	out, cert := certify(3, false, true)
	if cert != nil {
		t.Fatal("unexpected certificate")
	}
	s := string(encodeJSON(out, false))
	for _, want := range []string{`"polynomial":[1,2]`, `"lower_size":72`, `"interval_size":28`,
		`"rank_vector":[1,8,12,6,1]`, `"terminal":true`, `"mu":0`} {
		if !strings.Contains(s, want) {
			t.Errorf("summary %s lacks %s", s, want)
		}
	}
}

func TestAffineD4Certificate(t *testing.T) {
	_, cert := certify(4, true, false)
	s := string(cert)
	if !strings.HasPrefix(s, `{"format_version":1,"group":{"rank":5,"edges":[[0,1],[0,2],[0,3],[0,4]]}`) {
		t.Fatalf("unexpected header %.80s", s)
	}
	if !strings.Contains(s, `"polynomial":[1,3,2],"mu":2`) || !strings.HasSuffix(s, "}\n") {
		t.Fatal("summary of the D4 certificate is wrong")
	}
}

// TestCPythonSetOrder pins the set-table emulation (insertClean, resize and
// the merge paths of pyset.go) independently of the committed certificate:
// the expected order is CPython 3.14.7's list(g.lower(w)) for the original
// targeted.py; the 20 elements cross a table resize.
func TestCPythonSetOrder(t *testing.T) {
	g := star(2)
	w := g.elt([]int{1, 2, 0, 1, 2})
	want := [][]int{{0}, {0, 1}, {2, 1}, {1}, {1, 0, 2, 0, 1}, {2, 0, 1}, {0, 2}, {2, 1, 0}, {0, 2, 1}, {1, 0, 2, 1},
		{2, 0}, {1, 0, 2, 0}, {2}, {0, 1, 0}, {1, 0, 2}, {1, 0}, {0, 2, 0}, {}, {2, 0, 1, 0}, {0, 2, 0, 1}}
	got := g.lower(w).order()
	if len(got) != len(want) {
		t.Fatalf("lower(w) has %d elements, want %d", len(got), len(want))
	}
	for i, z := range got {
		if !slices.Equal(g.word(z), want[i]) {
			t.Fatalf("element %d of lower(w) has word %v, want %v", i, g.word(z), want[i])
		}
	}
}
