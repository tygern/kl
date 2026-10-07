# Leading Kazhdan–Lusztig coefficients in type E

A research project extending the methods of [Tyson Gern's type-D thesis](https://arxiv.org/abs/1304.6074) to ordinary, equal-parameter Kazhdan–Lusztig polynomials with a fully commutative lower endpoint.

The manuscript establishes:

- **Finite exceptional types:** the leading coefficient is 0 or 1 in E6, E7 and E8 for every upper endpoint. Together with earlier A and D results, this gives the bound for every finite simply laced Coxeter group.
- **Arbitrarily large ranks:** explicit full-support terminal reflection families in E(4r+1), r ≥ 3, with vanishing non-cover leading coefficients and eventual complete vanishing at each fixed rank.
- **Further examples:** complete leading-coefficient vanishing for an affine E8 family, and a complementary E10 construction.

The finite classification is computer-assisted; the uniform construction has a symbolic proof. Historical priority remains qualified. These results concern leading coefficients, not formulas for every coefficient of the polynomials or arbitrary upper endpoints in infinite type. Combinatorial-invariance transfer is a separate conditional observation.

Read the [manuscript PDF](output/pdf/exceptional-leading.pdf) or [LaTeX source](results/exceptional-leading.tex). The [v0.1.0 release](https://github.com/tygern/kl/releases/tag/v0.1.0) supplies a compact, self-contained proof supplement and checksums.

## Reproduce

Download and extract `exceptional-leading-proof.zip` from the release, enter `exceptional-leading-proof`, and run:

```sh
python3 run_proofs.py
```

Python 3.10+ and an installed C++17 compiler suffice; no network or Python packages are needed. Add `--full` to regenerate the historical full-group E6/E7 searches. Both modes passed after isolated extraction. The supplement includes all required inputs, expected outputs, provenance and a SHA-256 manifest; its README states exactly what is rebuilt and what is verified from saved certificates.

The repository retains computation sources and research reports in `computations/`, `research/` and `results/`. [Referee reports](research/journal-review/) document Astra-assisted mathematical, literature and editorial reviews; they are not journal decisions. Build the supplement with `python3 supplement/build_package.py`; rebuild the PDFs with `python3 research/build_pdfs.py` using an installed TeX toolchain.
