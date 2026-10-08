package d6

import (
	"go/ast"
	"go/parser"
	"go/token"
	"strings"
	"testing"
)

func perm(t *testing.T, w ...int64) Perm {
	p, err := NewPerm(w)
	if err != nil {
		t.Fatal(err)
	}
	return p
}

func TestPolyArithmetic(t *testing.T) {
	if !Equal(Add(nil, []int64{1, 1}, 1, 2), []int64{0, 2, 2}) {
		t.Fatal("add")
	}
	if !Equal(Mul([]int64{1, 1}, []int64{1, 1}), []int64{1, 2, 1}) {
		t.Fatal("mul")
	}
	if len(Trim([]int64{0, 0})) != 0 {
		t.Fatal("trim")
	}
}

func TestA2(t *testing.T) {
	g, _ := NewSparseCoxeter("A", 2)
	e := perm(t, 1, 2, 3)
	w0 := perm(t, 3, 2, 1)
	if g.Length(w0) != 3 || len(g.Lower(w0)) != 6 || !g.Leq(e, w0) {
		t.Fatal("A2 basics")
	}
	if !Equal(g.RPoly(e, w0), []int64{-1, 2, -2, 1}) {
		t.Fatalf("R(e,w0)=%v", g.RPoly(e, w0))
	}
}

func TestD6Root(t *testing.T) {
	g, _ := NewSparseCoxeter("D", 6)
	x := perm(t, -1, -2, 4, 3, 6, 5)
	w := perm(t, -1, -6, 3, -4, 5, -2)
	if !g.Leq(x, w) {
		t.Fatal("x <= w expected")
	}
	if len(g.Descents(w)) != 4 {
		t.Fatalf("descents %v", g.Descents(w))
	}
}

// The package must not contain a KL evaluator.
func TestNoEvaluator(t *testing.T) {
	fset := token.NewFileSet()
	pkgs, err := parser.ParseDir(fset, ".", nil, 0)
	if err != nil {
		t.Fatal(err)
	}
	banned := map[string]bool{"kl": true, "withdescent": true, "corrections": true,
		"verifyreciprocity": true, "klvalues": true, "mu": true}
	for _, pkg := range pkgs {
		for fname, f := range pkg.Files {
			if strings.HasSuffix(fname, "_test.go") {
				continue
			}
			ast.Inspect(f, func(n ast.Node) bool {
				if fd, ok := n.(*ast.FuncDecl); ok && banned[strings.ToLower(fd.Name.Name)] {
					t.Errorf("evaluator function %s in %s", fd.Name.Name, fname)
				}
				return true
			})
		}
	}
}
