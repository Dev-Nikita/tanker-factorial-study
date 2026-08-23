# Missing and ambiguous inputs

Nothing in this list has been invented. Each item states what is missing, what the
code does instead, and what has to be supplied before the corresponding number can
be called publication-ready.

Source codes: **SRC-A** = original tanker manuscript (.docx); **SRC-B** = 3³ factorial
brief (.docx); **SRC-C** = `red_hnups,+13.pdf`; **SRC-D** = *Теорія автомобіля* lecture notes.

---

> **Status update after the Mathcad archive became available.** Items M-02, M-03 and
> M-04 are now **resolved**: the archived ETAL worksheet supplies every parameter and
> every equation directly, so the manuscript model contains no calibrated quantity. See
> `ARCHIVE_TRACEABILITY.md`. The entries below are kept for the record, with their
> current status marked. M-01, M-05 and M-06 concern the separate 3^3 factorial study.

---

## M-01 — Physical meaning of factor x3 (factorial study only)

SRC-B defines x3 as "швидкість рушання" at 4.0 / 4.5 / 5.0 m s^-1. A vehicle starting
from rest has zero speed, so the term is ambiguous. Confirmed by the author as the
**target start-off speed**, which is how it is implemented. This item affects the
separate factorial study only; the manuscript prescribes its own reference manoeuvre.

---

## M-02 — RESOLVED: all parameters recovered from the archive

| ID | Parameter | Former status | Now |
|----|-----------|---------------|-----|
| P-1 | tractor mass | identified | **not applicable** - the archived model is a single three-body chain (chassis 5900 + body 700 + cargo 5000 kg); there is no separate tractor |
| P-2 | wheel radius `Rk` | identified as 0.5376 m | **archive: 0.530 m** (region 2102), independently confirmed by the vehicle data sheet |
| P-3 | driveline inertia `Itr` | identified as 3596 kg m^2 | **archive: 2025 kg m^2** (region 2112) |
| P-4 | tangential wheel damping `D_kch` | added term, identified | **removed** - the archive has no such term and none is needed; reproduction is within 1.14 % without it |

The earlier identification is superseded. Its only lasting value is that the identified
wheel radius, 0.5376 m, came within 1.4 % of the archived 0.530 m.

---

## M-03 — RESOLVED: the chain is chassis -> body -> cargo

SRC-A Table 1 reports three forces (`F_kch`, `F_chk`, `F_chv`). The archived worksheet
confirms the chain: `Mch` -> `Mkuz` -> `Mv`, coupled by `CCHK` and `CKV`. The two-tank
reading used in an earlier draft was wrong and has been corrected throughout.

---

## M-04 — RESOLVED: the source formulation is planar, and the archive supersedes the DOCX

An earlier version of this file asserted that the display object for eq. (1) renders as
"a different and much larger spatial model". That was an over-reading of the symbol set
and is withdrawn: the formulation is planar in the vertical longitudinal plane, with
`X`, `Z` and `phi` the longitudinal displacement, vertical displacement and pitch.

The question is now moot for the manuscript, because the equations are taken from the
archived worksheet rather than from the DOCX rendering. The implemented model retains
the longitudinal freedoms only and is described as a *reduced-order longitudinal model
derived from* the archived planar formulation.

---

## M-05 — Test-rig masses (factorial study only)

SRC-B gives tank geometry (50 L, ⌀370 mm, 500 mm) but not the dry tank mass, the
roller-carriage mass, the spring rates of the rig, or the liquid.

**Implemented:** liquid = water (user-confirmed); rig masses carry the label
`PLACEHOLDER_NOT_FOR_PUBLICATION` in `config/tanks.yaml`. The rig-scale factorial is
therefore reported as a *design study for the planned physical campaign*, never as a
result about the rig.

---

## M-06 — No physical measurements exist

No sensor data have been recorded. `data/experimental/` is empty and the validation
stage reports `NOT_AVAILABLE`. The manuscript states explicitly that the 27-run
factorial is a **computational** experiment. The CSV schema and the run protocol for
the future physical campaign are in `data/experimental/README.md` and
`docs/EXPERIMENT_PROTOCOL.md`.

---

## M-07 — Allowable longitudinal travel of the movable units

Needed to state whether an optimised configuration is feasible. Currently a
placeholder (0.25 m) in `config/springs.yaml`; used only to *report* whether the
travel limit is exceeded, never to constrain a published result.

---

## M-08 — Uncertainty magnitudes

`config/uncertainty.yaml` uses engineering tolerance classes (spring-rate class,
damper force tolerance, guide-condition range from SRC-A). These are justified but
not measured; the Monte-Carlo results are labelled accordingly.

---

## M-09 — NEW: the published first-gear peak forces cannot be reproduced

SRC-A Table 1 gives peak longitudinal forces for start-off *in first gear* at cargo
masses 1250–5000 kg. The archive-locked model exceeds them by 48–58 %. The archive holds
no first-gear traction characteristic — the ETAL worksheet contains a single 19-point
table whose peak is 65.5 kN — so the scenario cannot be reconstructed.

**Implemented:** those values are excluded from every comparison, the discrepancy is
reported in the manuscript as a limitation, and no parameter was adjusted to approach
them.

**Required:** the worksheet used for the first-gear runs, if it exists.
