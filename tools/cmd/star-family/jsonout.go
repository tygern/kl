package main

import (
	"bytes"
	"fmt"
	"strconv"
	"strings"
)

// A minimal ordered JSON writer reproducing Python's json.dumps byte for byte:
// indent=2 (every array element on its own line) or compact separators.

type kv struct {
	K string
	V any
}

type object []kv

func ints(v []int) []any {
	out := make([]any, len(v))
	for i, x := range v {
		out[i] = x
	}
	return out
}

func ints64(v []int64) []any {
	out := make([]any, len(v))
	for i, x := range v {
		out[i] = x
	}
	return out
}

func encodeJSON(v any, pretty bool) []byte {
	var b bytes.Buffer
	enc(&b, v, pretty, 0)
	return b.Bytes()
}

func enc(b *bytes.Buffer, v any, pretty bool, level int) {
	nl := func(l int) {
		if pretty {
			b.WriteByte('\n')
			b.WriteString(strings.Repeat("  ", l))
		}
	}
	switch t := v.(type) {
	case int:
		b.WriteString(strconv.Itoa(t))
	case int64:
		b.WriteString(strconv.FormatInt(t, 10))
	case bool:
		if t {
			b.WriteString("true")
		} else {
			b.WriteString("false")
		}
	case string:
		b.WriteString(strconv.Quote(t)) // ASCII text without special characters
	case []any:
		if len(t) == 0 {
			b.WriteString("[]")
			return
		}
		b.WriteByte('[')
		for i, x := range t {
			if i > 0 {
				b.WriteByte(',')
			}
			nl(level + 1)
			enc(b, x, pretty, level+1)
		}
		nl(level)
		b.WriteByte(']')
	case object:
		if len(t) == 0 {
			b.WriteString("{}")
			return
		}
		b.WriteByte('{')
		for i, e := range t {
			if i > 0 {
				b.WriteByte(',')
			}
			nl(level + 1)
			b.WriteString(strconv.Quote(e.K))
			b.WriteByte(':')
			if pretty {
				b.WriteByte(' ')
			}
			enc(b, e.V, pretty, level+1)
		}
		nl(level)
		b.WriteByte('}')
	default:
		panic(fmt.Sprintf("unsupported JSON value %T", v))
	}
}
