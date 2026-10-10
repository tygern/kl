# Leading Kazhdan–Lusztig coefficients with fully commutative lower endpoints in type E

Tyson Gern (Initial Capacity, tyson@initialcapacity.io, [ORCID 0009-0003-8288-8786](https://orcid.org/0009-0003-8288-8786))

## Results

For every fully commutative $x$ and arbitrary $w$ in $E_6$, $E_7$ and $E_8$, $\mu(x,w)\in\{0,1\}$. Green proved the type-$A$ case, and [Gern proved the type-$D$ case](https://arxiv.org/abs/1304.6074). The bound therefore holds for every finite simply laced Coxeter group, as Green asked.

The finite proof uses an exhaustive classification of terminal elements. Beyond Gern's type-$D$ list and its parabolic embeddings, there are two new elements: $w_7^E$ of length 28 and $w_8^E$ of length 50. By the maximum-descent lemma, each terminal has one eligible FC lower endpoint. For the new cases, the length gap is even and $\mu=0$; the only remaining KL value is the known $D_6$ coefficient.

Explicit infinite families of full-support terminals in affine $E_8$, $E_{10}$ and $E_{4r+1}$ for $r\geq3$ satisfy $\mu(x,b)=0$ for every FC $x$. The proofs hold for all parameter values, with finite computational inputs stated in Appendix A. For the affine family, five explicit braid witnesses persist in reduced words under a length-additive translation factorization. The FC catalogue and all-base-cover computation remain supplementary checks. The bound fails in affine $D_4$, where $\mu=2$, and without the FC hypothesis, where $\mu=10$ already in $E_6$.

## Read the paper

- [LaTeX source](results/exceptional-leading.tex) and matching [PDF](output/pdf/exceptional-leading.pdf) for release v0.6.1.
- [Citation audit](results/citation-audit.md): source locators, corrections, retained references, and verification limits.
- Section 2 develops star reduction, the maximum-descent lemma, and the geometric terminal test.
- Section 3 proves the finite theorem, then describes the supplementary structure of the terminals.
- Section 4 proves the reflection-cover lemmas, motivates the root constructions, and works through an $E_{13}$ example.
- Section 5 contains complementary computations and open problems.
- Appendix A explains the computational dependencies, proves the catalogue criterion, and prints the real-root witnesses.

Start with the finite proof, which uses the inherited maximum-descent argument and the new exceptional terminal classification by exhaustive parabolic pruning. The supplementary KL tables are optional for that proof. Then study the infinite families.

The [proof-study companion](results/proof-study-companion.md) contains small examples, eight exercises, and solution sketches on star transport, the Fan–Green cell argument, parabolic pruning, and the root construction.

## Reproduce the current revision

From the repository root:

```sh
go run -C supplement ./cmd/package
```

Extract `supplement/dist/exceptional-leading-proof.zip`, enter `exceptional-leading-proof`, and run:

```sh
go run ./cmd/proofs -full
```

For the finite theorem alone, run `go run ./cmd/proofs -finite`. The runner regenerates the classifications, verifies the printed terminal data, and recomputes the residual $D_6$ coefficient. Its summary contains only those inputs. The `-finite` and `-full` flags are mutually exclusive.

Go 1.22 or later suffices. The module uses only the standard library. The runner checks the manifest, rebuilds the computations in an isolated working tree, compares regenerated outputs, and checks the mathematical summary. Without `-full`, it omits the heavier cross-checks and compares against saved full-group $E_6/E_7$ snapshots. See the [supplement README](supplement/README.md) for the precise scope and [development guide](supplement/DEVELOPING.md) for individual commands.

[Release v0.6.1](https://github.com/tygern/kl/releases/tag/v0.6.1) contains the revised prose, two Dynkin diagrams, and the [repository-wide style audit](results/style-audit-2026-10-10.json), with a matching PDF and proof supplement. [Release v0.6.0](https://github.com/tygern/kl/releases/tag/v0.6.0) preserves the preceding structural proof revision. [Release v0.5.1](https://github.com/tygern/kl/releases/tag/v0.5.1) preserves the preceding citation revision. Each archive records its manuscript and source hashes.

The [repository tools](tools/README.md) reproduce the additional research certificates and build the PDF:

```sh
go run -C tools ./cmd/research-checks
go run -C tools ./cmd/build-pdfs
go run -C tools ./cmd/check-pdf
```

The PDF commands require an installed TeX toolchain and Poppler. Computational certificates are in `results/` and the corresponding `research/` directories. Historical notes under `research/` record earlier stages and may contain superseded statements; the manuscript and current supplement define the present claims.

## AI assistance and licensing

OpenAI and Anthropic models assisted with research, proof drafts, code, and internal checks. The author directed the work and takes responsibility for the manuscript. Internal AI reviews are not peer review; [PROVENANCE.md](PROVENANCE.md) records the workflow and the distinction between computational checks and review claims.

Code is licensed under the [MIT License](LICENSE). The manuscript and PDFs are licensed under [CC BY 4.0](LICENSE-TEXT.md).
