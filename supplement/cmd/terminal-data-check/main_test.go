package main

import (
	"os"
	"path/filepath"
	"testing"
)

func TestPoincareE6Order(t *testing.T) {
	var sum int64
	for _, c := range poincare([]int{2, 5, 6, 8, 9, 12}) {
		sum += c
	}
	if sum != 51840 {
		t.Fatalf("|W(E6)| = %d, want 51840", sum)
	}
}

func TestRootBounds(t *testing.T) {
	edges := [][2]int{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {2, 5}}
	n, m, err := rootBounds(6, edges)
	if err != nil || n != 72 || m != 3 {
		t.Fatalf("E6 roots: %d %d %v", n, m, err)
	}
	if got := independentSets(6, edges); got != 22 {
		t.Fatalf("independent sets %d", got)
	}
}

// dataRoot returns the directory holding results/e6-independent-certificate.json:
// the repository root (supplement/cmd/<x> -> ../../..), or payload/ inside the
// shipped archive (<archive>/cmd/<x> -> ../../payload; ../../../payload is also tried).
func dataRoot(t *testing.T) string {
	t.Helper()
	for _, c := range []string{"../../..", "../../payload", "../../../payload"} {
		if _, err := os.Stat(filepath.Join(c, "results", "e6-independent-certificate.json")); err == nil {
			return c
		}
	}
	t.Skip("certificates not found at ../../../results, ../../payload/results or ../../../payload/results")
	return ""
}

func TestRepositoryData(t *testing.T) {
	out, err := run(dataRoot(t))
	if err != nil {
		t.Fatal(err)
	}
	if out.Status != "All assertions passed" {
		t.Fatal(out.Status)
	}
}
