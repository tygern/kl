# Leading Kazhdan–Lusztig coefficients with fully commutative lower endpoints in type E

Tyson Gern (Initial Capacity, tyson@initialcapacity.io, [ORCID 0009-0003-8288-8786](https://orcid.org/0009-0003-8288-8786))

## Results

The paper proves that for every fully commutative $x$ and arbitrary $w$ in the Weyl groups of types $E_6$, $E_7$ and $E_8$ the leading coefficient $\mu(x,w)$ of the ordinary Kazhdan–Lusztig polynomial is $0$ or $1$. Combined with Green's theorem for type $A$ (proved for affine $A$ as well) and [Gern's type-$D$ thesis](https://arxiv.org/abs/1304.6074), this settles Green's 2009 question (Chmutov's "Green's 0–1 conjecture") for every finite simply laced Coxeter group. The proof reduces, by Green's star-operation argument and Fan's monomial cell theory, to a finite list of noncommuting terminal elements, which is computed exhaustively: it is Gern's type-$D$ list plus exactly two new elements, $w_7$ in $E_7$ (length 28) and $w_8$ in $E_8$ (length 50), and every new element has even gap to its unique eligible lower endpoint, so no Kazhdan–Lusztig computation beyond Gern's $D_6$ value is needed. The paper also constructs infinite families of full-support terminal elements in $E_{4r+1}$ ($r\ge 3$), in affine $E_8$ (lengths $27+92k$) and in $E_{10}$ for which $\mu(x,b)=0$ for every fully commutative $x$ and every member of the family, and shows that the bound is sharp: it fails in affine $D_4$ ($\mu=2$) and it fails without the fully commutative hypothesis ($\mu=10$ already in $E_6$). To our knowledge the exceptional bound, the terminal lists and the infinite families have not appeared before.

The finite classification is computer-assisted; the families have symbolic proofs. A complete $E_6$ table and the affine $E_8$ odd-gap coefficients of a second terminal family are reported in the paper's final section (Section 5), and the structure of $w_7$ and $w_8$ in Remark 3.3; all are certified in the supplement.

## Use of AI tools

This work was carried out with substantial assistance from large language model systems. OpenAI's GPT-6.1 Sol and GPT-6 Astra (the "Astra agents" of the internal notes), run as teams of agents, performed the literature searches, proposed and drafted the constructions and proofs, wrote the enumeration and verification code, produced internal review notes, and drafted the text. Anthropic's Claude Fable 5.1 then carried out an independent adversarial review that re-derived every lemma against the primary sources, recomputed every computational claim with separately written code, and found the strengthenings and corrections incorporated in the current version. The Go port of the proof supplement (release v0.3.0) was written by Anthropic's Claude Sonnet 5.5 and verified by Claude Fable 5.1 against the earlier programs' outputs. The author directed the work, checked the proofs, computations and references, and takes full responsibility for the content. The internal review notes are machine-generated and are not peer review. The complete workflow, including the prompts in `prompts/`, is described in [PROVENANCE.md](PROVENANCE.md).

## Licensing

Code is licensed under the [MIT License](LICENSE). The manuscript text and the PDFs are licensed under [CC BY 4.0](LICENSE-TEXT.md).

## Contents

- Paper: [LaTeX source](results/exceptional-leading.tex) and [PDF](output/pdf/exceptional-leading.pdf).
- Independent review (7 October 2026, Claude Fable 5.1): [report](research/claude-review-2026-10-07/REVIEW.md), [findings with verifier verdicts](research/claude-review-2026-10-07/findings-digest.md); the reviewers' code is in the repository history at commit d997526.
- Archived research notes (superseded, folded into Section 5 of the paper): [research/archive/notes/](research/archive/notes/).
- Internal AI review and research notes: [research/ai-review-notes/](research/ai-review-notes/), [research/editorial/](research/editorial/) and the other `research/` directories. Every text note carries a banner stating that it is machine-generated and is not peer review or a journal decision.
- Proof supplement: the Go module [supplement/](supplement/DEVELOPING.md) (one command per program, the proof runner and the archive packager); the archive README and validation records under `supplement/` (built archive under `output/supplement/`, not tracked).
- Certificates and data: `results/` and the `research/` directories.
- Repository tools: the Go module [tools/](tools/README.md) (PDF build and checks, and the research-note certificates that are not part of the proof supplement, with the runner `research-checks`).
- The repository contains Go only. The Python and C++ programs of releases v0.1.0 and v0.2.0 are in release v0.2.0 and in the repository history at commit d997526.

## Reproduce

[Release **v0.4.0**](https://github.com/tygern/kl/releases/tag/v0.4.0) of the repository (archive `exceptional-leading-proof.zip`, `SHA256SUMS`, release notes) certifies every computational claim of the current paper with a single Go module. The earlier releases remain the historical tags: [v0.3.0](https://github.com/tygern/kl/releases/tag/v0.3.0) is the first Go-module supplement (v0.4.0 adds six certified steps and removes the last Python and C++ files from the repository), [v0.2.0](https://github.com/tygern/kl/releases/tag/v0.2.0) shipped the same computations as Python and C++ programs (the Go port was verified against their outputs certificate for certificate), and [v0.1.0](https://github.com/tygern/kl/releases/tag/v0.1.0) is the earlier version of the manuscript. Extract the supplement archive, enter `exceptional-leading-proof`, and run

```sh
go run ./cmd/proofs
```

Go 1.22 or later suffices; the module uses the standard library alone, so no network access, module download, package or compiler other than Go is needed. Add `-full` to regenerate the full-group $E_6/E_7$ searches and the heavier certificates. The supplement includes all inputs, expected outputs, provenance and a SHA-256 manifest; its README states exactly what is rebuilt and what is verified from saved certificates.

Within the repository, `go run -C supplement ./cmd/package` builds the supplement archive, `go build -C supplement -o bin/ ./cmd/...` builds every command into `supplement/bin/` (run from the repository root, a command rewrites its committed certificate in place, or prints it for the commands that write to stdout; see `supplement/DEVELOPING.md` for the redirection form), and `go run -C tools ./cmd/build-pdfs` followed by `go run -C tools ./cmd/check-pdf` rebuilds and checks the PDF with an installed TeX toolchain and poppler (see `tools/README.md`).
