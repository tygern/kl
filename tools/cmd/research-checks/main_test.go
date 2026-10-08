package main

import (
	"strings"
	"testing"
)

func parse(t *testing.T, text string) any {
	t.Helper()
	v, err := decodeJSON([]byte(text))
	if err != nil {
		t.Fatal(err)
	}
	return v
}

func TestDiffJSON(t *testing.T) {
	cases := []struct{ a, b, want string }{
		{`{"x": 1, "y": [1, 2.0]}`, `{"y": [1.0, 2], "x": 1}`, ""},
		{`{"x": 1}`, `{"x": 2}`, "$.x: 1 vs 2"},
		{`{"x": [1, 2]}`, `{"x": [1]}`, "$.x: list lengths 2 vs 1"},
		{`{"x": 1}`, `{"x": 1, "z": true}`, "$.z: missing in the committed file"},
		{`{"s": "a"}`, `{"s": "b"}`, `$.s: "a" vs "b"`},
		{`[true, null]`, `[true, null]`, ""},
	}
	for _, c := range cases {
		if got := diffJSON(parse(t, c.a), parse(t, c.b), "$"); got != c.want {
			t.Errorf("diff(%s, %s) = %q, want %q", c.a, c.b, got, c.want)
		}
	}
}

func TestTable(t *testing.T) {
	out := table([]row{{"i6-shape", "results/i6-shape-certificate.json", "equal", "", 0.5}, {"x", "y", "differs", "$.k: 1 vs 2", 1}})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 4 || !strings.HasPrefix(lines[0], "step") || !strings.Contains(lines[1], "equal") || strings.TrimSpace(lines[3]) != "$.k: 1 vs 2" {
		t.Fatalf("unexpected table:\n%s", out)
	}
}
