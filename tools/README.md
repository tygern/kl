# Repository tools

The Go module `github.com/tygern/kl/tools` builds and checks the manuscript PDF and regenerates the research-note certificates. The proof supplement's programs are in `supplement/`. Go 1.22 or later suffices; the module uses only the standard library.

Run commands from the repository root with `go run -C tools ./cmd/<name>`. Commands that read or write repository files search parent directories for `results/` and `tools/go.mod` to locate the repository root (`-root DIR` overrides this). `go build -C tools -o bin/ ./cmd/...` writes binaries to the gitignored `tools/bin/` directory; binaries are not tracked.

| Command | What it does | Output |
|---|---|---|
| `research-checks` | Runs the six commands below in a scratch tree under `tmp/research-checks/`, compares every regenerated file with the committed one (JSON semantically, text byte for byte; no key is ignored), prints a table and exits non-zero on any difference. It generates the E9 FC catalogue with the supplement's `fc-catalogue -rank 9` first, since `catalogue-eligible` reads it. | table on stdout |
| `sources-manifest` | SHA-256 digests and source URLs of the locally archived PDFs under `sources/` (a gitignored local archive; nothing is downloaded). | `sources/manifest.json` |
| `i6-shape` | Shape of Gern's D6 interval `[x, w]` with `x = (-1,-2,4,3,6,5)`, `w = (-1,-6,3,-4,5,-2)`: 1,676 elements, rank 11, 12 atoms (so not a principal lower interval) and not a lattice; a second part re-derives the ideal by direct subwords with no order comparison. | `results/i6-shape-certificate.json`, `results/i6-independent-subword.json` |
| `exceptional-cosets` | Root certificates for the E7 and E8 cosets containing copies of Gern's D6 interval: lengths, stripped right descents and the images of the two remaining adjacent simple roots. | text on stdout, committed as `research/exceptional-cosets.txt` |
| `star-family` | Kazhdan-Lusztig polynomial of the star family (a central generator joined to `m` commuting leaves) for `m = 1..8`, `P_{x,w} = sum_k (C(m,k) - C(m,k-1)) q^k`, with the lower-ideal and interval sizes `4^m + 2^m` and `3^m + 1`; for `m = 4` (affine D4) every KL dependency is saved. | `research/broad_cells/star_family_checks.json`, `research/broad_cells/affine_d4_kl_certificate.json` |
| `affine-d-checks` | The lower Bruhat interval of `(0,1,3,4,2)^2` in affine D4 (rank vector and the four rank-3 edges obstructing a principal finite-D unfolding), and the lower covers of Gern's `w` in D_n, `n = 4, 6, ..., 24`. | `research/broad_affine/unfolding_obstruction.json`, `research/broad_affine/gern_coatoms.json` |
| `catalogue-eligible -rank N` | For each full-support terminal `w` of the E_N conjugate search, the fully commutative `x <= w` with at least the descents of `w` (the eligible lower endpoints). Reads the complete FC catalogue `research/en_independent/eN-fc.json`, a gitignored output of `go run -C supplement ./cmd/fc-catalogue -rank N`. | `research/en_families/eN_full_support_eligible.json` (committed for `N = 9`) |
| `review-finite` | Independent review checks of every printed finite terminal word, descents, support independence and orthogonal-reflection descriptors. | Assertions on stdout |
| `review-infinite` | Independent review checks of affine inversion classes, braid witnesses, root-lowering sequences and exact-integer samples of the uniform identities. For proofs at arbitrary parameter values, see the manuscript. | Assertions on stdout |
| `build-pdfs` | Builds `results/*.tex` with `latexmk` into `tmp/pdfs/build/<stem>/` and copies the PDF to `output/pdf/`. | `output/pdf/<stem>.pdf`, `results/pdf-build.log`, `results/pdf-build-manifest.json` |
| `check-pdf` | Checks every PDF under `output/pdf/` (page count, title and author, text on every page, a 65 pt margin test on the word boxes) and writes page contact sheets. | `results/pdf-qa.json`, `tmp/pdfs/<stem>-contact-N.png` |

`build-pdfs` and `check-pdf` need `latexmk` (TeX Live) and poppler (`pdfinfo`, `pdftotext`, `pdftoppm`) on `PATH`; everything else needs only Go.

```sh
go run -C tools ./cmd/research-checks      # a few seconds; rewrites nothing in place
go run -C tools ./cmd/build-pdfs && go run -C tools ./cmd/check-pdf
cd tools && gofmt -l . && go vet ./... && go test ./...
```

To rewrite a committed certificate after a deliberate change, run its command directly (for example `go run -C tools ./cmd/star-family`). Then run `research-checks` to compare the committed file with the program's output.
