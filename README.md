# primes-strata-gue-isomorphism

**The Music of the Primes Meets the Rhythm of the Strata: A Same-Universality-Class Test Between Riemann Zeta Zeros and Real Milankovitch Cyclostratigraphy**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)

**Author:** Ruqing Chen · GUT Geoservice Inc., Montreal · ruqing@hotmail.com

A **companion** to Part I of a fifth racetrack applying random matrix theory to Earth's
spatial rhythms: *Spatial Level Repulsion in the Stratigraphic Record*
([Zenodo 20774581](https://zenodo.org/records/20774581)).

---

## The Question

The nontrivial zeros of the Riemann zeta function follow the **GUE** (Gaussian Unitary
Ensemble) of random matrix theory — the Montgomery–Odlyzko law. Part I showed that the
down-core spacings of Milankovitch cycles in real cyclostratigraphic sections also reach
the GUE class. Placed on the same scale, with **purely real data on both sides**, are
they statistically isomorphic?

## Result

| | Riemann zeros | Gubbio | Blue Lias | GUE theory |
|---|---|---|---|---|
| spacings *n* | 4519 | 31 | 98 | — |
| ⟨r⟩ | 0.617 | 0.735 | 0.748 | 0.603 |
| Δ⟨r⟩ vs zeros | — | +0.12 | +0.13 | — |
| KS *p* vs zeros | — | 0.65 | **0.031** | — |

- The 4519 real zero spacings track GUE theory almost exactly (⟨r⟩ = 0.617 ≈ 0.603).
- **Both** real strata sections (the two of Part I) are systematically **more rigid**
  than the zeros (Δ⟨r⟩ = +0.12 and +0.13) — closer to a perfect lattice, exactly as
  astronomical quasi-periodic forcing predicts.
- Gubbio alone (31 spacings) cannot resolve this (KS p = 0.65, small sample). The larger
  **Blue Lias** section (98 spacings) **rejects** equality with the zeros (p = 0.031), as
  do the two combined (n = 129, p = 0.037).
- That two unrelated sections — different basins, ages, lithologies — show the **same**
  rigidity offset indicates it is a general property of astronomically forced strata.

![Q-Q and CDF](figures/qq_isomorphism.pdf)

## Honest Conclusion

Both belong to the **GUE universality class** — a real and deep statistical-physics
correspondence, the same kind Montgomery & Dyson found between primes and atomic nuclei.
But they are **not the same distribution**: two independent strata sections are
systematically more rigid, significantly so once enough spacings are included. Sharing a
universality class does **not** imply any causal or ontological link between primes and
strata; universality is precisely that microscopically unrelated systems can share
macroscopic statistics.

> **Same key (GUE), not the same melody.** To sharpen the comparison, the next step is a
> longer real stratigraphic record, not more zeros.

This is a companion to, **not** Part II of, the racetrack. Part II is reserved for the
structural fault-spacing analog, which tests the spatial-shadow hypothesis directly.

## Data (real, traceable)

- `data/riemann_zeros_real.csv` — 4520 real Riemann zeros (imaginary parts γₙ, t ∈ (14, 5000]),
  computed by **David Platt** to 100-bit precision, archived in the **LMFDB**
  ([lmfdb.org/zeros/zeta](https://www.lmfdb.org/zeros/zeta/), a database of 103.8 billion zeros).
  Decoded from the native Platt binary format; first five values match the published
  Odlyzko table (14.134725, 21.022040, 25.010858, 30.424876, 32.935062).
- `data/gubbio_real.csv` — real Gubbio Contessa magnetic susceptibility
  (Sinnesael et al. 2016, [PANGAEA 864450](https://doi.org/10.1594/PANGAEA.864450), CC-BY-3.0).
- `data/bluelias_full.csv` — real Blue Lias magnetic susceptibility + lithology
  (Weedon et al. 2019, [PANGAEA 896875](https://doi.org/10.1594/PANGAEA.896875), CC-BY-4.0).

Both are plain CSV (openable in Excel, R, pandas, or any text editor).

## Reproduce

```bash
pip install -r requirements.txt
cd code
python run_isomorphism.py                 # the same-universality-class test
python make_figure.py                     # regenerate the Q-Q + CDF figure
python decode_platt_zeros.py <file.dat>   # decode a raw Platt/LMFDB binary file
```

## The Platt Decoder

`code/decode_platt_zeros.py` decodes David Platt's native binary zero format
(40-byte header + 13-byte incremental fixed-point records in units of 2⁻¹⁰¹). It is a
reusable tool for anyone working with the LMFDB's 103.8-billion-zero database, with
built-in verification against the known Odlyzko values and a monotonicity check.

## Files

```
primes-strata-gue-isomorphism/
├── README.md / LICENSE / requirements.txt / CITATION.cff / .zenodo.json
├── code/
│   ├── run_isomorphism.py     one-click reproduction (reads CSV)
│   ├── make_figure.py         regenerate the figure
│   └── decode_platt_zeros.py  Platt/LMFDB binary decoder
├── data/
│   ├── riemann_zeros_real.csv 4520 real zeros (Platt/LMFDB)
│   ├── gubbio_real.csv        real Gubbio (PANGAEA 864450)
│   └── bluelias_full.csv      real Blue Lias (PANGAEA 896875)
├── figures/qq_isomorphism.pdf
└── paper/
    ├── paper.pdf              6-page paper
    ├── paper.tex
    └── figs/
```

## Related

- **Part I (published):** Spatial Level Repulsion in the Stratigraphic Record. Zenodo: https://zenodo.org/records/20774581
- RMT program (Racetracks 1–4): zenodo 20766310, 20768130, 20768420, 20768849

## References

- Montgomery, H.L. (1973). The pair correlation of zeros of the zeta function.
- Odlyzko, A.M. (1987). On the distribution of spacings between zeros of the zeta function.
- Platt, D. & Trudgian, T. (2021). The Riemann hypothesis is true up to 3·10¹².
- The LMFDB Collaboration (2024). Zeros of ζ(s). https://www.lmfdb.org/zeros/zeta/

## License

MIT (code). Real data redistributed under their sources' terms (LMFDB; PANGAEA CC-BY-3.0).
