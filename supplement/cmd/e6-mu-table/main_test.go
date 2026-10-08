package main

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// TestMain lets the test binary act as the command when E6MU_RUN is set.
func TestMain(m *testing.M) {
	if args := os.Getenv("E6MU_RUN"); args != "" {
		os.Args = append([]string{"e6-mu-table"}, strings.Fields(args)...)
		main()
		os.Exit(0)
	}
	os.Exit(m.Run())
}

// TestA4 runs the full pipeline on A4 (|W| = 120): all assertions inside the
// command (left/right recursion agreement, nonnegativity, degree bound) must
// pass, and the certificate must carry the known group data.
func TestA4(t *testing.T) {
	out := filepath.Join(t.TempDir(), "a4.json")
	exe, err := os.Executable()
	if err != nil {
		t.Fatal(err)
	}
	cmd := exec.Command(exe)
	cmd.Env = append(os.Environ(), "E6MU_RUN=-type A -rank 4 -out "+out)
	if b, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("run failed: %v\n%s", err, b)
	}
	b, err := os.ReadFile(out)
	if err != nil {
		t.Fatal(err)
	}
	s := string(b)
	for _, want := range []string{`"group_order": 120`, `"longest_length": 10`, `"positive_roots": 10`, `"status": "passed"`} {
		if !strings.Contains(s, want) {
			t.Errorf("certificate lacks %s", want)
		}
	}
}
