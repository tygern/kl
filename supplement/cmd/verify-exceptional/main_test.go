package main

import (
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

func TestPoincare(t *testing.T) {
	if sum(polynomial([]int{1, 4, 5, 7, 8, 11})) != 51840 {
		t.Fatal("E6 order")
	}
	if sum(polynomial([]int{1, 5, 7, 9, 11, 13, 17})) != 2903040 {
		t.Fatal("E7 order")
	}
}

func TestSigned(t *testing.T) {
	got := signed([]int64{1})
	if got[0] != -2 || got[1] != -1 {
		t.Fatalf("unexpected %v", got)
	}
}

// readJSON exits the process on failure, so the rejection case runs in a
// re-executed copy of the test binary.
func TestReadJSONRejectsTrailingData(t *testing.T) {
	if os.Getenv("VE_TRAILING_CHILD") == "1" {
		readJSON("bad.json")
		os.Exit(0)
	}
	dir := t.TempDir()
	if err := os.WriteFile(filepath.Join(dir, "ok.json"), []byte("{\"a\": 1}\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(dir, "bad.json"), []byte("{\"a\": 1}}"), 0o644); err != nil {
		t.Fatal(err)
	}
	cmd := exec.Command(os.Args[0], "-test.run=TestReadJSONRejectsTrailingData")
	cmd.Dir = dir
	cmd.Env = append(os.Environ(), "VE_TRAILING_CHILD=1")
	if err := cmd.Run(); err == nil {
		t.Fatal("trailing data after JSON value was accepted")
	}
}

func TestReadJSONAcceptsTrailingWhitespace(t *testing.T) {
	dir := t.TempDir()
	old, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	if err := os.Chdir(dir); err != nil {
		t.Fatal(err)
	}
	defer os.Chdir(old)
	if err := os.WriteFile("ok.json", []byte("{\"a\": 1}\n  \n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if m := readJSON("ok.json"); m["a"] == nil {
		t.Fatal("valid JSON not read")
	}
}
