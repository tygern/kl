package main

import (
	"os"
	"path/filepath"
	"testing"
)

func TestPoincare(t *testing.T) {
	// Order of E8 = product of degrees; Poincare polynomial at 1.
	p := poincare([]int{2, 8, 12, 14, 18, 20, 24, 30})
	var sum int64
	for _, c := range p {
		sum += c
	}
	if sum != 696729600 {
		t.Fatalf("|W(E8)| = %d", sum)
	}
	if len(p) != 121 {
		t.Fatalf("length %d", len(p))
	}
}

func writeCert(t *testing.T, body string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "cert.json")
	if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func panics(f func()) (r bool) {
	defer func() {
		if recover() != nil {
			r = true
		}
	}()
	f()
	return false
}

// The longest element of the A2 parabolic on {0,1} has both generators as
// descents but is not weak-terminal: w(alpha_0+alpha_1) is negative.
func TestWeakOrderConditionRejected(t *testing.T) {
	p := writeCert(t, `{"bad":[{"word":[0,1,0],"length":3,"Rmask":3,"Lmask":3,"bottoms":[]}]}`)
	if !panics(func() { checkCertificate(p, 6) }) {
		t.Fatal("non weak-terminal element accepted")
	}
}

// A single generator adjacent to nothing it dominates is weak-terminal.
func TestWeakTerminalAccepted(t *testing.T) {
	p := writeCert(t, `{"bad":[{"word":[0],"length":1,"Rmask":1,"Lmask":1,"bottoms":[]}]}`)
	if panics(func() { checkCertificate(p, 6) }) {
		t.Fatal("valid terminal word rejected")
	}
}
