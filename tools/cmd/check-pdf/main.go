// Command check-pdf reads the built PDFs under output/pdf/, checks them and
// writes page contact sheets and results/pdf-qa.json.
//
// It ports research/check_pdf.py (Python with Pillow, pypdf and pdfplumber)
// to Go with the standard library only. External programs used (poppler):
//
//	pdfinfo     page count, Title and Author
//	pdftotext   -bbox: word boxes for the text and margin checks
//	pdftoppm    -png -r 60: page renders in tmp/pdfs/<stem>-<page>.png
//
// Contact sheets are built with image, image/draw and image/png.
//
// Margin test. The Python original tested every CHARACTER box from
// pdfplumber (x0 < 65 or x1 > page width - 65). pdftotext -bbox only reports
// WORD boxes, so this port tests word boxes (xMin < 65 or xMax > width - 65).
// A word box is the union of its character boxes, so any character outside
// the margins puts its word outside too: the word test flags every page the
// character test flags (it can only be stricter, never weaker, apart from
// differences between pdftotext and pdfplumber glyph-advance conventions).
// A page "has text" when pdftotext finds at least one non-empty word on it.
//
// Run from the repository root (go run -C tools ./cmd/check-pdf runs in tools/,
// so the repository root is found by walking up to the directory that contains
// tools/go.mod; -root overrides). Imports no other package of this module.
package main

import (
	"bytes"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"image"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"unicode/utf8"
)

const (
	marginPt     = 65.0
	renderDPI    = "60"
	visualNote   = "Contact sheets and representative full page inspected separately."
	sheetTiles   = 6
	tilesPerRow  = 3
	tileW, tileH = 430, 560
)

// report is one entry of results/pdf-qa.json (key order as the original).
type report struct {
	PDF               string   `json:"pdf"`
	Pages             int      `json:"pages"`
	Bytes             int64    `json:"bytes"`
	Title             *string  `json:"title"`
	Author            *string  `json:"author"`
	AllPagesHaveText  bool     `json:"all_pages_have_text"`
	WithinPageMargins bool     `json:"within_page_margins"`
	ContactSheets     []string `json:"contact_sheets"`
	VisualInspection  string   `json:"visual_inspection"`
}

func main() {
	rootFlag := flag.String("root", "", "repository root (default: nearest ancestor of the working directory containing tools/go.mod)")
	flag.Parse()
	if err := run(*rootFlag); err != nil {
		fmt.Fprintln(os.Stderr, "check-pdf:", err)
		os.Exit(1)
	}
}

func findRoot(flagValue string) (string, error) {
	start := flagValue
	if start == "" {
		wd, err := os.Getwd()
		if err != nil {
			return "", err
		}
		start = wd
	}
	start, err := filepath.Abs(start)
	if err != nil {
		return "", err
	}
	if resolved, err := filepath.EvalSymlinks(start); err == nil {
		start = resolved
	}
	if flagValue != "" {
		if !isRoot(start) {
			return "", fmt.Errorf("%s is not a repository root (needs results/ and tools/go.mod)", flagValue)
		}
		return start, nil
	}
	for dir := start; ; {
		if isRoot(dir) {
			return dir, nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", errors.New("repository root not found (no ancestor contains results/ and tools/go.mod); use -root")
		}
		dir = parent
	}
}

// isRoot reports whether dir is the repository root: it holds results/ and
// the tools module file tools/go.mod.
func isRoot(dir string) bool {
	if st, err := os.Stat(filepath.Join(dir, "results")); err != nil || !st.IsDir() {
		return false
	}
	st, err := os.Stat(filepath.Join(dir, "tools", "go.mod"))
	return err == nil && st.Mode().IsRegular()
}

func run(rootFlag string) error {
	root, err := findRoot(rootFlag)
	if err != nil {
		return err
	}
	for _, tool := range []string{"pdfinfo", "pdftotext", "pdftoppm"} {
		if _, err := exec.LookPath(tool); err != nil {
			return fmt.Errorf("required program %s not found on PATH (install poppler)", tool)
		}
	}
	pdfs, err := filepath.Glob(filepath.Join(root, "output", "pdf", "*.pdf"))
	if err != nil {
		return err
	}
	sort.Strings(pdfs)
	reports := []report{}
	for _, pdf := range pdfs {
		r, err := checkOne(root, pdf)
		if err != nil {
			return err
		}
		reports = append(reports, r)
	}
	out, err := encodeReports(reports)
	if err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Join(root, "results"), 0o755); err != nil {
		return err
	}
	if err := os.WriteFile(filepath.Join(root, "results", "pdf-qa.json"), out, 0o644); err != nil {
		return err
	}
	_, err = os.Stdout.Write(out)
	return err
}

// encodeReports mimics Python json.dumps(indent=2, ensure_ascii=True) + "\n".
func encodeReports(reports []report) ([]byte, error) {
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	enc.SetIndent("", "  ")
	if err := enc.Encode(reports); err != nil {
		return nil, err
	}
	return asciiEscape(buf.Bytes()), nil
}

// asciiEscape rewrites every non-ASCII rune as \uXXXX (surrogate pairs above
// the BMP), as Python's json module does by default.
func asciiEscape(in []byte) []byte {
	var out bytes.Buffer
	for len(in) > 0 {
		r, size := utf8.DecodeRune(in)
		in = in[size:]
		switch {
		case r < 0x80:
			out.WriteByte(byte(r))
		case r < 0x10000:
			fmt.Fprintf(&out, "\\u%04x", r)
		default:
			v := r - 0x10000
			fmt.Fprintf(&out, "\\u%04x\\u%04x", 0xd800+(v>>10), 0xdc00+(v&0x3ff))
		}
	}
	return out.Bytes()
}

func rel(root, path string) (string, error) {
	r, err := filepath.Rel(root, path)
	if err != nil {
		return "", err
	}
	return filepath.ToSlash(r), nil
}

func checkOne(root, pdf string) (report, error) {
	stem := strings.TrimSuffix(filepath.Base(pdf), ".pdf")
	relPDF, err := rel(root, pdf)
	if err != nil {
		return report{}, err
	}
	info, err := readInfo(pdf)
	if err != nil {
		return report{}, err
	}
	if info.pages <= 0 {
		return report{}, fmt.Errorf("%s: no pages", relPDF)
	}
	pages, err := readWordBoxes(pdf)
	if err != nil {
		return report{}, err
	}
	if len(pages) != info.pages {
		return report{}, fmt.Errorf("%s: pdfinfo reports %d pages but pdftotext produced %d", relPDF, info.pages, len(pages))
	}
	var outside []string
	for number, page := range pages {
		if len(page.words) == 0 {
			return report{}, fmt.Errorf("%s: page %d has no text", relPDF, number+1)
		}
		// The manuscript is set with 0.9in (64.8pt) margins; text may
		// overhang the text block by a few points.
		for _, w := range page.words {
			if w.xMin < marginPt || w.xMax > page.width-marginPt {
				outside = append(outside, fmt.Sprintf("page %d %q x0=%v x1=%v", number+1, w.text, w.xMin, w.xMax))
			}
		}
	}
	if len(outside) > 0 {
		return report{}, fmt.Errorf("%s: words outside the %gpt margins:\n  %s", relPDF, marginPt, strings.Join(outside, "\n  "))
	}

	tmp := filepath.Join(root, "tmp", "pdfs")
	if err := os.MkdirAll(tmp, 0o755); err != nil {
		return report{}, err
	}
	if err := removeStale(tmp, stem); err != nil {
		return report{}, err
	}
	cmd := exec.Command("pdftoppm", "-png", "-r", renderDPI, pdf, filepath.Join(tmp, stem))
	if msg, err := cmd.CombinedOutput(); err != nil {
		return report{}, fmt.Errorf("pdftoppm failed on %s: %v\n%s", relPDF, err, msg)
	}
	renders, err := numberedRenders(tmp, stem)
	if err != nil {
		return report{}, err
	}
	if len(renders) != info.pages {
		return report{}, fmt.Errorf("%s: %d page renders for %d pages", relPDF, len(renders), info.pages)
	}
	thumbs := make([]*image.RGBA, 0, len(renders))
	for i, path := range renders {
		im, err := loadRGBA(path)
		if err != nil {
			return report{}, err
		}
		thumbs = append(thumbs, makeTile(im, i+1))
	}
	contacts := []string{}
	for start := 0; start < len(thumbs); start += sheetTiles {
		end := start + sheetTiles
		if end > len(thumbs) {
			end = len(thumbs)
		}
		sheet := makeSheet(thumbs[start:end])
		target := filepath.Join(tmp, fmt.Sprintf("%s-contact-%d.png", stem, start/sheetTiles+1))
		if err := savePNG(target, sheet); err != nil {
			return report{}, err
		}
		r, err := rel(root, target)
		if err != nil {
			return report{}, err
		}
		contacts = append(contacts, r)
	}
	st, err := os.Stat(pdf)
	if err != nil {
		return report{}, err
	}
	return report{
		PDF: relPDF, Pages: info.pages, Bytes: st.Size(),
		Title: info.title, Author: info.author,
		AllPagesHaveText: true, WithinPageMargins: true,
		ContactSheets: contacts, VisualInspection: visualNote,
	}, nil
}

var renderName = regexp.MustCompile(`^(.*)-(\d+)\.png$`)

// removeStale deletes earlier numbered renders and contact sheets of the stem
// so that leftovers from a longer document cannot be mistaken for pages.
func removeStale(dir, stem string) error {
	entries, err := os.ReadDir(dir)
	if err != nil {
		return err
	}
	contact := regexp.MustCompile("^" + regexp.QuoteMeta(stem) + `-contact-\d+\.png$`)
	for _, e := range entries {
		name := e.Name()
		m := renderName.FindStringSubmatch(name)
		if contact.MatchString(name) || (m != nil && m[1] == stem) {
			if err := os.Remove(filepath.Join(dir, name)); err != nil {
				return err
			}
		}
	}
	return nil
}

// numberedRenders returns <stem>-<n>.png paths sorted by the integer n.
func numberedRenders(dir, stem string) ([]string, error) {
	entries, err := os.ReadDir(dir)
	if err != nil {
		return nil, err
	}
	type numbered struct {
		n    int
		path string
	}
	var found []numbered
	for _, e := range entries {
		m := renderName.FindStringSubmatch(e.Name())
		if m == nil || m[1] != stem {
			continue
		}
		n, err := strconv.Atoi(m[2])
		if err != nil {
			return nil, err
		}
		found = append(found, numbered{n, filepath.Join(dir, e.Name())})
	}
	sort.Slice(found, func(i, j int) bool { return found[i].n < found[j].n })
	paths := make([]string, len(found))
	for i, f := range found {
		paths[i] = f.path
	}
	return paths, nil
}
