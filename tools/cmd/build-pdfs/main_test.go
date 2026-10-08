package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestManifestFormat(t *testing.T) {
	got, err := marshalManifest([]record{{
		Source: "results/a.tex", PDF: "output/pdf/a.pdf", SourceSHA256: "aa", PDFSHA256: "bb",
		Bytes: 3, Command: []string{"latexmk", "-pdf"}, CWD: "repository root",
	}})
	if err != nil {
		t.Fatal(err)
	}
	want := `[
  {
    "source": "results/a.tex",
    "pdf": "output/pdf/a.pdf",
    "source_sha256": "aa",
    "pdf_sha256": "bb",
    "bytes": 3,
    "command": [
      "latexmk",
      "-pdf"
    ],
    "cwd": "repository root"
  }
]
`
	if string(got) != want {
		t.Fatalf("manifest mismatch:\n%s", got)
	}
}

func TestNormalizeNewlines(t *testing.T) {
	if got := string(normalizeNewlines([]byte("a\r\nb\rc\n"))); got != "a\nb\nc\n" {
		t.Fatalf("got %q", got)
	}
}

func TestFindRootAndRel(t *testing.T) {
	dir := t.TempDir()
	if err := os.MkdirAll(filepath.Join(dir, "results"), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.MkdirAll(filepath.Join(dir, "tools", "cmd"), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(dir, "tools", "go.mod"), []byte("module x\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	old, _ := os.Getwd()
	defer os.Chdir(old)
	if err := os.Chdir(filepath.Join(dir, "tools", "cmd")); err != nil {
		t.Fatal(err)
	}
	root, err := findRoot("")
	if err != nil {
		t.Fatal(err)
	}
	r, err := rel(root, filepath.Join(root, "tmp", "pdfs", "build", "x"))
	if err != nil || r != "tmp/pdfs/build/x" {
		t.Fatalf("rel = %q, %v", r, err)
	}
	if _, err := rel(root, filepath.Dir(root)); err == nil || !strings.Contains(err.Error(), "outside") {
		t.Fatalf("expected outside-repository error, got %v", err)
	}
}
