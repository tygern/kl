package main

import (
	"bytes"
	"encoding/json"
	"fmt"
	"strconv"
	"strings"
)

// kv / obj form an insertion-ordered JSON object, so that the certificate
// keeps the key order of the original program.
type kv struct {
	K string
	V any
}
type obj []kv

func jsonString(s string) string {
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	if err := enc.Encode(s); err != nil {
		panic(err)
	}
	return strings.TrimSuffix(buf.String(), "\n")
}

// render writes v as JSON in the style of Python's json.dumps: with indent < 0
// the compact one-line form (", " and ": " separators), otherwise one space
// of indentation per level.
func render(sb *strings.Builder, v any, indent, level int) {
	nl := func(l int) {
		if indent >= 0 {
			sb.WriteString("\n")
			sb.WriteString(strings.Repeat(" ", indent*l))
		}
	}
	sep := ", "
	if indent >= 0 {
		sep = ","
	}
	list := func(n int, item func(i int)) {
		if n == 0 {
			sb.WriteString("[]")
			return
		}
		sb.WriteString("[")
		for i := 0; i < n; i++ {
			if i > 0 {
				sb.WriteString(sep)
			}
			nl(level + 1)
			item(i)
		}
		nl(level)
		sb.WriteString("]")
	}
	switch x := v.(type) {
	case nil:
		sb.WriteString("null")
	case bool:
		if x {
			sb.WriteString("true")
		} else {
			sb.WriteString("false")
		}
	case int:
		sb.WriteString(strconv.Itoa(x))
	case string:
		sb.WriteString(jsonString(x))
	case []int:
		list(len(x), func(i int) { sb.WriteString(strconv.Itoa(x[i])) })
	case []string:
		list(len(x), func(i int) { sb.WriteString(jsonString(x[i])) })
	case [][]int:
		list(len(x), func(i int) { render(sb, x[i], indent, level+1) })
	case []any:
		list(len(x), func(i int) { render(sb, x[i], indent, level+1) })
	case obj:
		if len(x) == 0 {
			sb.WriteString("{}")
			return
		}
		sb.WriteString("{")
		for i, e := range x {
			if i > 0 {
				sb.WriteString(sep)
			}
			nl(level + 1)
			sb.WriteString(jsonString(e.K))
			sb.WriteString(": ")
			render(sb, e.V, indent, level+1)
		}
		nl(level)
		sb.WriteString("}")
	default:
		panic(fmt.Sprintf("render: unsupported type %T", v))
	}
}

func toJSON(v any, indent int) string {
	var sb strings.Builder
	render(&sb, v, indent, 0)
	return sb.String()
}
