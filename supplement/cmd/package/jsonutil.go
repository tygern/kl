package main

// Small JSON helpers: insertion-ordered objects (encoding/json sorts map
// keys), exact number handling (json.Number, never float64) and a semantic
// comparison that reports the first difference with its JSON path. The
// runner (cmd/proofs) carries an identical copy; the two commands share no
// code, as the Python originals did not.

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

// object is a JSON object that marshals its keys in insertion order.
type object struct {
	keys   []string
	values map[string]any
}

func newObject() *object { return &object{values: map[string]any{}} }

func (o *object) set(key string, value any) *object {
	if _, ok := o.values[key]; !ok {
		o.keys = append(o.keys, key)
	}
	o.values[key] = value
	return o
}

func (o *object) MarshalJSON() ([]byte, error) {
	var b bytes.Buffer
	b.WriteByte('{')
	for i, k := range o.keys {
		if i > 0 {
			b.WriteByte(',')
		}
		kb, err := json.Marshal(k)
		if err != nil {
			return nil, err
		}
		b.Write(kb)
		b.WriteByte(':')
		vb, err := json.Marshal(o.values[k])
		if err != nil {
			return nil, err
		}
		b.Write(vb)
	}
	b.WriteByte('}')
	return b.Bytes(), nil
}

// decodeJSON parses a complete JSON document, keeping numbers exact.
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

// encodeJSON renders a value with two-space indentation and a final newline,
// the layout of Python's json.dumps(value, indent=2) + '\n'.
func encodeJSON(v any) ([]byte, error) {
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(v); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

// normalize re-decodes a value so that only map[string]any, []any,
// json.Number, string, bool and nil remain (ordered objects become maps,
// Go integers become json.Number).
func normalize(v any) (any, error) {
	data, err := json.Marshal(v)
	if err != nil {
		return nil, err
	}
	return decodeJSON(data)
}

func numberEqual(a, b json.Number) bool {
	if a == b {
		return true
	}
	x, okx := new(big.Rat).SetString(string(a))
	y, oky := new(big.Rat).SetString(string(b))
	return okx && oky && x.Cmp(y) == 0
}

// diffJSON compares two normalized values and returns "" when they are
// semantically equal (numbers by value, objects regardless of key order,
// lists in order), otherwise a description of the first difference.
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
				return p + ": missing on the first side"
			}
			if !iny {
				return p + ": missing on the second side"
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

// sameJSON reports whether two (possibly unnormalized) values are
// semantically equal; the second result is the first difference.
func sameJSON(a, b any) (bool, string) {
	na, err := normalize(a)
	if err != nil {
		return false, err.Error()
	}
	nb, err := normalize(b)
	if err != nil {
		return false, err.Error()
	}
	d := diffJSON(na, nb, "$")
	return d == "", d
}

// Accessors for decoded JSON. Each one panics with a failure describing the
// missing or mistyped field, the counterpart of a Python KeyError/TypeError.

func field(v any, key string) any {
	m, ok := v.(map[string]any)
	if !ok {
		panic(failure{fmt.Sprintf("expected a JSON object with key %q, found %T", key, v)})
	}
	x, ok := m[key]
	if !ok {
		panic(failure{fmt.Sprintf("missing key %q", key)})
	}
	return x
}

func items(v any) []any {
	l, ok := v.([]any)
	if !ok {
		panic(failure{fmt.Sprintf("expected a JSON list, found %T", v)})
	}
	return l
}

func mapping(v any) map[string]any {
	m, ok := v.(map[string]any)
	if !ok {
		panic(failure{fmt.Sprintf("expected a JSON object, found %T", v)})
	}
	return m
}

func sortedKeys(m map[string]any) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

func integer(v any) int64 {
	n, ok := v.(json.Number)
	if !ok {
		panic(failure{fmt.Sprintf("expected an integer, found %T", v)})
	}
	x, err := n.Int64()
	if err != nil {
		panic(failure{fmt.Sprintf("expected an integer, found %s", n)})
	}
	return x
}

func boolean(v any) bool {
	b, ok := v.(bool)
	if !ok {
		panic(failure{fmt.Sprintf("expected a boolean, found %T", v)})
	}
	return b
}

// pick copies the named keys of an object, in the given order.
func pick(v any, keys ...string) *object {
	o := newObject()
	for _, k := range keys {
		o.set(k, field(v, k))
	}
	return o
}

// literal decodes a JSON literal written in the source.
func literal(text string) any {
	v, err := decodeJSON([]byte(text))
	if err != nil {
		panic(failure{"bad literal " + text + ": " + err.Error()})
	}
	return v
}

// failure is raised (by panic) for any failed assertion; main recovers it.
type failure struct{ msg string }

func (f failure) Error() string { return f.msg }
