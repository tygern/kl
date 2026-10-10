# Leading Kazhdan–Lusztig coefficients with fully commutative lower endpoints in type E

Tyson Gern (Initial Capacity, tyson@initialcapacity.io, [ORCID 0009-0003-8288-8786](https://orcid.org/0009-0003-8288-8786))

## Results

For every fully commutative $x$ and arbitrary $w$ in $E_6$, $E_7$ and $E_8$, the paper proves $\mu(x,w)\in\{0,1\}$. Together with Green's type-$A$ theorem and [Gern's type-$D$ thesis](https://arxiv.org/abs/1304.6074), this settles Green's question for every finite simply laced Coxeter group.

The finite proof reduces to an exhaustive classification of terminal elements. Beyond Gern's type-$D$ list and its parabolic embeddings, there are two new elements: $w_7^E$ of length 28 and $w_8^E$ of length 50. A maximum-descent lemma forces one eligible FC lower endpoint for each terminal. Parity eliminates the new cases; the only remaining KL value is the known $D_6$ coefficient.

The paper also constructs infinite families of full-support terminals in affine $E_8$, $E_{10}$ and $E_{4r+1}$ for $r\geq3$, with $\mu(x,b)=0$ for every FC $x$. These arguments are symbolic in their parameters, with finite computational inputs stated in Appendix A. The affine family propagates five explicit braid witnesses through a length-additive translation factorization. The FC catalogue and all-base-cover computation remain supplementary checks. The bound fails in affine $D_4$, where $\mu=2$, and without the FC hypothesis, where $\mu=10$ already in $E_6$.

## Read the paper

- [LaTeX source](results/exceptional-leading.tex) and [PDF](output/pdf/exceptional-leading.pdf), synchronized for release v0.6.0.
- [Citation audit](results/citation-audit.md): source locators, corrections, retained references, and verification limits.
- Section 2 develops star reduction, the maximum-descent lemma, and the geometric terminal test.
- Section 3 proves the finite theorem, then gives the supplementary structure of the terminals.
- Section 4 proves the reflection-cover lemmas, motivates the root constructions, and works through an $E_{13}$ example.
- Section 5 gives complementary computations and open problems.
- Appendix A explains the computational dependencies, proves the catalogue criterion, and prints the real-root witnesses.

For a first reading, follow the finite proof before studying the infinite families. The key inputs are the inherited maximum-descent argument and the new exceptional terminal classification by exhaustive parabolic pruning; the large supplementary KL tables are not prerequisites.

The [proof-study companion](results/proof-study-companion.md) works through star transport, the Fan–Green cell argument, parabolic pruning, and the root construction with small examples, eight exercises, and solution sketches.

## Reproduce the current revision

From the repository root:

```sh
go run -C supplement ./cmd/package
```

Extract `supplement/dist/exceptional-leading-proof.zip`, enter `exceptional-leading-proof`, and run:

```sh
go run ./cmd/proofs -full
```

For just the finite theorem, run `go run ./cmd/proofs -finite` instead. This regenerates the classifications, verifies the printed terminal data, and recomputes the residual $D_6$ coefficient, with a summary restricted to those inputs. The `-finite` and `-full` flags are mutually exclusive.

Go 1.22 or later suffices. The module uses only the standard library. The runner checks the manifest, rebuilds the computations in an isolated working tree, compares regenerated outputs, and checks the mathematical summary. Without `-full`, it omits the heavier cross-checks and compares against saved full-group $E_6/E_7$ snapshots. See the [supplement README](supplement/README.md) for the precise scope and [development guide](supplement/DEVELOPING.md) for individual commands.

[Release v0.6.0](https://github.com/tygern/kl/releases/tag/v0.6.0) contains the affine braid-witness proof, compact descriptions of the exceptional terminals, and a finite-only verification mode, with the updated PDF, companion, audit, and supplement. [Release v0.5.1](https://github.com/tygern/kl/releases/tag/v0.5.1) preserves the preceding citation revision. Each archive records its manuscript and source hashes.

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
