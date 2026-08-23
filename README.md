# Longitudinal dynamics of a tank vehicle with segmented cargo mounting

Reproducible computational study supporting the manuscript *"Separating structural-load
mitigation from traction demand in tank vehicles with segmented cargo mounting"*.

The model is built **directly from the archived Mathcad worksheet** that produced the
original published results. Every parameter, function and the tabulated traction
characteristic are extracted programmatically; **no parameter is calibrated**.

---

## Quick start

```bash
pip install -r requirements.txt
python3 scripts/run_all.py
```

Runtime is under a minute. The pipeline extracts the archive, runs studies A1–A7,
regenerates the eight manuscript figures and runs the 28 automated tests.

Outputs:

- `results/archive/` — extracted scalars, functions and file hashes
- `results/archive_study/` — result tables and archived time histories
- `results/archive_figures/` — the eight manuscript figures (vector PDF + 600 dpi PNG)

---

## Scope

This repository serves two separate studies.

**Manuscript pipeline** — `tools/parse_mathcad_archive.py`, `scripts/11`, `scripts/12`,
`src/tanker_dynamics/archive_model.py`, `archive_studies.py`. This is the study described
in `manuscript/`.

**A separate 3³ full-factorial study** — `scripts/02`–`10`, `docs/EXPERIMENT_PROTOCOL.md`
and the older `model.py`/`solver.py`/`inverse.py`. It targets a future *experimental*
paper on a scaled rig and is deliberately excluded from the manuscript pipeline. Its code
is retained unchanged.

---

## What the study does

1. **Extracts the archive.** The ETAL worksheet yields 72 scalar definitions, 63 function
   definitions and a 19-point traction table, each tagged with its worksheet region id.
   See `docs/ARCHIVE_TRACEABILITY.md`.
2. **Reproduces the published results** with no fitted parameter.
3. **Compares rigid and segmented mounting twice**: forward, with the archived traction
   characteristic (the published, uncontrolled comparison), and inverse, with the same
   prescribed start-off imposed on both (the controlled one).
4. **Sweeps the parameter space** — cargo and body mass, coupling stiffness, guide
   friction, damping ratio, manoeuvre duration and three reference trajectories.

### Headline results

**Reproduction, zero calibration.** Maximum error over nine published quantities: 1.14 %.

| Quantity | Published | Model | Error |
|---|---:|---:|---:|
| Terminal speed, m/s | 13.058 | 13.0576 | −0.003 % |
| Distance, m | 7582 | 7582.8 | +0.011 % |
| Engine work, MJ | 53.37 | 52.76 | −1.14 % |
| Braking time, s | 25.054 | 24.967 | −0.35 % |
| Braking distance, m | 163.575 | 163.59 | +0.011 % |
| Braking peak force, N | 6800 | 6828 | +0.41 % |

**Fixed powertrain (forward).** Peak wheel force falls 7.7–20.5 % with segmented
mounting; largest at the highest cargo mass. The published claim reproduces in direction.

**Matched manoeuvre (inverse).** Nothing falls:

| Quantity | Rigid | Segmented | Change |
|---|---:|---:|---:|
| Peak traction force, kN | 17.70 | 17.71 | −0.04 % |
| Peak traction power, kW | 63.25 | 63.30 | −0.09 % |
| Traction work, MJ | 0.6912 | 0.6913 | −0.01 % |
| Minimum effective mass, t | 11.60 | 10.40 | +10.4 % |

Over the whole examined domain the peak-force reduction is never positive, reaching
−27 % near `T/T_n ≈ 2`. The traction work is larger in every case, as the work–energy
balance requires for a passive coupling.

---

## Layout

```
tools/parse_mathcad_archive.py   XMCD extraction (scalars, functions, SHA-256)
src/tanker_dynamics/
  archive_model.py               ETAL worksheet transcribed; region ids in the docstring
  archive_studies.py             forward integration and matched-manoeuvre inverse dynamics
  plotting.py                    figure style
  model.py, solver.py, doe.py    the separate factorial study
scripts/11, 12, run_all          manuscript pipeline
scripts/00-10                    the separate factorial study
docs/ARCHIVE_TRACEABILITY.md     every parameter and equation mapped to a worksheet region
docs/MISSING_INPUTS.md           what is still unresolved
tests/                           28 tests
```

---

## Verification

`python3 -m pytest tests -q` — 28 tests. The ones that matter most:

- **archive fidelity** — the scalars in `ArchiveParams` equal the extracted worksheet
  values, and the derived quantities satisfy `Mchk = 1/Mch + 1/Mkuz` etc.;
- **no calibrated parameters** — all optional damping terms default to zero;
- **momentum balance** — the three body equations sum to `F_kch − A1`, so the internal
  friction pairs cancel;
- **reproduction** — terminal speed, distance, braking time and braking distance against
  the published values;
- **stiff limit** — with `C → 5×10¹⁰` N/m the segmented model approaches the rigid one;
- **matched-manoeuvre work never decreases** — an initially unstrained passive coupling
  cannot lower the traction work;
- **damping extension** matches the formula declared in the manuscript, and all three
  reference trajectories satisfy their stated boundary conditions.

The last two tests each caught a real defect during development: a missing friction
reaction on the frame, and a sign error in the inverse-dynamics routine. Both had
produced an apparent benefit that does not exist.

Numerical convergence: tightening `rtol`/`atol` by an order of magnitude changes the
reported peaks by less than 0.05 %.

---

## Honesty constraints observed

- **No calibrated parameters.** Everything comes from the archived worksheet. The
  `D_kch` damping term used in an earlier draft is removed: it is not needed.
- **No mixing of worksheets.** The archive holds three different computational cases;
  only the ETAL worksheet is used. See `ARCHIVE_TRACEABILITY.md` §E.
- **No fabricated measurements.** The study is computational; nothing is simulated in
  place of data.
- **Failures reported.** The published first-gear peak forces are *not* reproduced
  (48–58 % discrepancy) and no parameter was adjusted to approach them.
- **Bounded claims.** The energy result is proved from the work–energy balance for
  initially unstrained passive couplings; the peak-force result is stated as empirical
  over the examined domain, not as a theorem.
- **Provenance published.** The manuscript cites the source study and the archive, and
  quotes the SHA-256 digest of the worksheet the model is built from.

---

## Licence

MIT. See `LICENSE`. If you use this work, cite the article; see `CITATION.cff`.
