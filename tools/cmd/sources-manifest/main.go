// Command sources-manifest records the locally archived public sources with
// SHA-256 digests.
//
// It ports research/source_manifest.py to Go (standard library only; no other
// package of this repository is imported). Nothing is downloaded: the program
// hashes the PDF files already present under sources/ and records the URLs
// they were obtained from, then writes sources/manifest.json with the same
// schema, key order and formatting as the Python original (two-space indent,
// trailing newline). sources/ is gitignored; only the PDFs' local presence is
// required.
//
// Invocation, from the repository root:
//
//	go run -C tools ./cmd/sources-manifest
//
// `go run -C tools` runs the program with tools/ as the working directory, so
// the repository root is found by walking upward from the working directory
// until a directory containing results/ and tools/go.mod is reached (the -root
// flag overrides this).
package main

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"os"
	"path/filepath"
)

// source is one archived file and the URL it was downloaded from; the order of
// the slice is the order of the original URLS dictionary.
type source struct {
	File string
	URL  string
}

var sources = []source{
	{"combinatorial-invariance-2026.pdf", "https://raw.githubusercontent.com/openai/math/main/preprints/Combinatorial-Invariance-of-Kazhdan-Lusztig-Polynomials-September-24-2026/paper.pdf"},
	{"gern-thesis-2013.pdf", "https://arxiv.org/pdf/1304.6074"},
	{"green-leading-2008.pdf", "https://arxiv.org/pdf/0801.1650"},
	{"green-jones-traces-2007.pdf", "https://arxiv.org/pdf/math/0509362"},
	{"green-losonczy-cells-2001.pdf", "https://arxiv.org/pdf/math/0102003"},
	{"jones-deodhar-2007.pdf", "https://arxiv.org/pdf/0711.1391"},
	{"woo-patterns-2006.pdf", "https://arxiv.org/pdf/math/0611328"},
}

// record is one entry of "sources"; field order is the key order of the original.
type record struct {
	File      string `json:"file"`
	SourceURL string `json:"source_url"`
	SizeBytes int    `json:"size_bytes"`
	SHA256    string `json:"sha256"`
}

type notArchived struct {
	Title  string `json:"title"`
	URL    string `json:"url"`
	Reason string `json:"reason"`
}

type manifest struct {
	ResearchDateLocal string        `json:"research_date_local"`
	Timezone          string        `json:"timezone"`
	Note              string        `json:"note"`
	Sources           []record      `json:"sources"`
	NotArchived       []notArchived `json:"not_archived"`
}

func main() {
	rootFlag := flag.String("root", "", "repository root (default: found by walking up from the working directory)")
	flag.Parse()
	if flag.NArg() != 0 {
		fatalf("unexpected arguments: %v", flag.Args())
	}
	root := *rootFlag
	if root == "" {
		var err error
		if root, err = findRoot(); err != nil {
			fatalf("%v", err)
		}
	}
	records := make([]record, 0, len(sources))
	for _, s := range sources {
		path := filepath.Join(root, "sources", s.File)
		data, err := os.ReadFile(path)
		if err != nil {
			fatalf("cannot read sources/%s: %v", s.File, err)
		}
		if !bytes.HasPrefix(data, []byte("%PDF-")) {
			fatalf("Not a PDF: sources/%s", s.File)
		}
		records = append(records, hashRecord(s, data))
	}
	out, err := marshalManifest(records)
	if err != nil {
		fatalf("%v", err)
	}
	if err := os.WriteFile(filepath.Join(root, "sources", "manifest.json"), out, 0o644); err != nil {
		fatalf("%v", err)
	}
	fmt.Printf("Archived source manifest: %d verified PDF files.\n", len(records))
}

func hashRecord(s source, data []byte) record {
	sum := sha256.Sum256(data)
	return record{
		File:      "sources/" + s.File,
		SourceURL: s.URL,
		SizeBytes: len(data),
		SHA256:    hex.EncodeToString(sum[:]),
	}
}

func marshalManifest(records []record) ([]byte, error) {
	m := manifest{
		ResearchDateLocal: "2026-10-06",
		Timezone:          "America/Chicago",
		Note:              "Locally archived snapshots, not a guarantee of external acceptance or later version stability.",
		Sources:           records,
		NotArchived: []notArchived{{
			Title:  "Chmutov 2014 thesis",
			URL:    "https://hdl.handle.net/2027.42/108802",
			Reason: "University PDF returned an access challenge; indexed primary text and author-uploaded text inspected; see literature audit.",
		}},
	}
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(m); err != nil { // Encode appends the trailing newline
		return nil, err
	}
	return buf.Bytes(), nil
}

// findRoot walks upward from the working directory to the first directory that
// contains both results/ and tools/go.mod.
func findRoot() (string, error) {
	dir, err := os.Getwd()
	if err != nil {
		return "", err
	}
	for {
		_, e1 := os.Stat(filepath.Join(dir, "results"))
		_, e2 := os.Stat(filepath.Join(dir, "tools", "go.mod"))
		if e1 == nil && e2 == nil {
			return dir, nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", errors.New("repository root not found (no directory with results/ and tools/go.mod above the working directory); use -root")
		}
		dir = parent
	}
}

func fatalf(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "sources-manifest: "+format+"\n", args...)
	os.Exit(1)
}
