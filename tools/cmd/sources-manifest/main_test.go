package main

import (
	"strings"
	"testing"
)

func TestHashRecord(t *testing.T) {
	r := hashRecord(source{"a.pdf", "https://example.org/a"}, []byte("%PDF-abc"))
	if r.File != "sources/a.pdf" || r.SizeBytes != 8 || len(r.SHA256) != 64 {
		t.Fatalf("unexpected record %+v", r)
	}
}

func TestManifestFormat(t *testing.T) {
	got, err := marshalManifest([]record{{File: "sources/a.pdf", SourceURL: "u&v", SizeBytes: 3, SHA256: "aa"}})
	if err != nil {
		t.Fatal(err)
	}
	s := string(got)
	for _, frag := range []string{
		"{\n  \"research_date_local\": \"2026-10-06\",\n  \"timezone\": \"America/Chicago\",\n",
		"  \"sources\": [\n    {\n      \"file\": \"sources/a.pdf\",\n      \"source_url\": \"u&v\",\n      \"size_bytes\": 3,\n      \"sha256\": \"aa\"\n    }\n  ],\n  \"not_archived\": [\n",
	} {
		if !strings.Contains(s, frag) {
			t.Fatalf("missing fragment %q in\n%s", frag, s)
		}
	}
	if !strings.HasSuffix(s, "}\n") {
		t.Fatalf("missing trailing newline")
	}
}

func TestSourceList(t *testing.T) {
	if len(sources) != 7 || sources[0].File != "combinatorial-invariance-2026.pdf" {
		t.Fatal("source list changed")
	}
}
