package main

// Semantic JSON comparison: numbers by value (json.Number, never float64),
// objects regardless of key order, lists in order. The first difference is
// reported with its JSON path.

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"math/big"
	"os"
	"sort"
	"strconv"
)

func decodeJSON(data []byte) (any, error) {
	dec := json.NewDecoder(bytes.NewReader(data))
	dec.UseNumber()
	var v any
	if err := dec.Decode(&v); err != nil {
		return nil, err
	}
	if _, err := dec.Token(); err != io.EOF {
		return nil, fmt.Errorf("trailing data after the JSON document")
	}
	return v, nil
}

func readJSON(path string) (any, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	v, err := decodeJSON(data)
	if err != nil {
		return nil, fmt.Errorf("%s: %v", path, err)
	}
	return v, nil
}

func numberEqual(a, b json.Number) bool {
	if a == b {
		return true
	}
	x, okx := new(big.Rat).SetString(string(a))
	y, oky := new(big.Rat).SetString(string(b))
	return okx && oky && x.Cmp(y) == 0
}

// diffJSON returns "" when a and b are semantically equal, otherwise a
// description of the first difference ("committed vs regenerated").
func diffJSON(a, b any, path string) string {
	switch x := a.(type) {
	case map[string]any:
		y, ok := b.(map[string]any)
		if !ok {
			return fmt.Sprintf("%s: object vs %T", path, b)
		}
		keys := make([]string, 0, len(x)+len(y))
		seen := map[string]bool{}
		for k := range x {
			keys = append(keys, k)
			seen[k] = true
		}
		for k := range y {
			if !seen[k] {
				keys = append(keys, k)
			}
		}
		sort.Strings(keys)
		for _, k := range keys {
			p := path + "." + k
			xv, inx := x[k]
			yv, iny := y[k]
			if !inx {
				return p + ": missing in the committed file"
			}
			if !iny {
				return p + ": missing in the regenerated file"
			}
			if d := diffJSON(xv, yv, p); d != "" {
				return d
			}
		}
		return ""
	case []any:
		y, ok := b.([]any)
		if !ok {
			return fmt.Sprintf("%s: list vs %T", path, b)
		}
		if len(x) != len(y) {
			return fmt.Sprintf("%s: list lengths %d vs %d", path, len(x), len(y))
		}
		for i := range x {
			if d := diffJSON(x[i], y[i], path+"["+strconv.Itoa(i)+"]"); d != "" {
				return d
			}
		}
		return ""
	case json.Number:
		y, ok := b.(json.Number)
		if !ok {
			return fmt.Sprintf("%s: number vs %T", path, b)
		}
		if !numberEqual(x, y) {
			return fmt.Sprintf("%s: %s vs %s", path, x, y)
		}
		return ""
	case string:
		y, ok := b.(string)
		if !ok {
			return fmt.Sprintf("%s: string vs %T", path, b)
		}
		if x != y {
			return fmt.Sprintf("%s: %q vs %q", path, x, y)
		}
		return ""
	case bool:
		y, ok := b.(bool)
		if !ok {
			return fmt.Sprintf("%s: bool vs %T", path, b)
		}
		if x != y {
			return fmt.Sprintf("%s: %v vs %v", path, x, y)
		}
		return ""
	case nil:
		if b != nil {
			return fmt.Sprintf("%s: null vs %T", path, b)
		}
		return ""
	}
	return fmt.Sprintf("%s: unsupported value of type %T", path, a)
}
