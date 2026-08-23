# Archive traceability

Every quantity in the manuscript model maps to a definition in the archived Mathcad
worksheet. Extraction is programmatic (`tools/parse_mathcad_archive.py`); the outputs
are `results/archive/archive_parameters.csv`, `archive_functions.csv` and
`archive_manifest.csv` (SHA-256 of every archived file).

**Authoritative source: `10 тест зад (РОЗГІН- РУХ-ГАЛЬМУВ) 23.10.09р. ETAL..xmcd`**
(Mathcad Professional 14.1). It is the only worksheet that covers the full
acceleration–motion–braking sequence described in the article.

## Provenance codes

| Code | Meaning |
|---|---|
| `ARCHIVE_ETAL` | value read from the ETAL worksheet |
| `ARCHIVE_VEHICLE` | value from `ПАРАМЕТРИ ZiL-131.doc` (cross-check only) |
| `DERIVED` | computed from `ARCHIVE_ETAL` values inside the worksheet |
| `NUMERICAL` | solver setting, no physical content |
| `CALIBRATED` | fitted — **none remain in the manuscript model** |

---

## A. Scalars

| ID | Symbol | Value | Unit | Region | Provenance |
|---|---|---|---|---|---|
| A-01 | `Rk` | 0.530 | m | 2102 | ARCHIVE_ETAL |
| A-02 | `g` | 9.81 | m/s² | 2107 | ARCHIVE_ETAL |
| A-03 | `Alf_d` | 0.06 | – | 2108 | ARCHIVE_ETAL |
| A-04 | `Mch` | 5900 | kg | 2109 | ARCHIVE_ETAL |
| A-05 | `Mkuz` | 700 | kg | 2110 | ARCHIVE_ETAL |
| A-06 | `Mv` | 5000 | kg | 2111 | ARCHIVE_ETAL |
| A-07 | `Itr` | 2025 | kg·m² | 2112 | ARCHIVE_ETAL |
| A-08 | `Alf_ch` | 0.02 | – | 2143 | ARCHIVE_ETAL |
| A-09 | `Alf_k` | 0.02 | – | 4065 | ARCHIVE_ETAL |
| A-10 | `CKCH` | 1.2×10⁶ | N·m/rad | 4069 | ARCHIVE_ETAL |
| A-11 | `CCHK` | 2.0×10⁵ | N/m | 4070 | ARCHIVE_ETAL |
| A-12 | `CKV` | 2.0×10⁵ | N/m | 4071 | ARCHIVE_ETAL |
| A-13 | `Mgal_max` | 1600 | N·m | 3803 | ARCHIVE_ETAL |
| A-14 | `T1_poch`, `T1_kin` | 0, 600 | s | 2675–2676 | ARCHIVE_ETAL |
| A-15 | `XCH1, Bk, LCHK` | 2.50, 2.35, 0.15 | m | 4072–4074 | ARCHIVE_ETAL |
| A-16 | `XK2, Bv, LKV` | 3.00, 2.85, 0.15 | m | 4075–4077 | ARCHIVE_ETAL |
| A-17 | `MSchkv` | 11600 | kg | 2140 | DERIVED = Mch+Mkuz+Mv |
| A-18 | `MSkv` | 5700 | kg | 2141 | DERIVED = Mkuz+Mv |
| A-19 | `Mchk` | 1/Mch+1/Mkuz | kg⁻¹ | 2113 | DERIVED |
| A-20 | `Mkv` | 1/Mkuz+1/Mv | kg⁻¹ | 2142 | DERIVED |
| A-21 | traction table `v1`,`pt1` | 19 points, max 65 472 N | m/s, N | 2103–2104 | ARCHIVE_ETAL |

Note on A-15/A-16: `XCH1 − Bk − LCHK = 0` and `XK2 − Bv − LKV = 0`, so the spring free
lengths cancel and the coupling forces reduce to `−CCHK·u1` and `−CKV·u2`, with `u`
measured from the unstretched state. This is a change of origin only.

---

## B. Functions

| ID | Worksheet | Region | Implementation |
|---|---|---|---|
| B-01 | `MprKCH(y) = CKCH*(y1 − y2/Rk)` | 4078 | `archive_model.rhs` |
| B-02 | `FprCHK(y) = CCHK*((XCH1 − y3) − Bk − LCHK)` | 4079 | `archive_model.rhs` |
| B-03 | `FprKV(y) = CKV*((XK2 − y4) − Bv − LKV)` | 4080 | `archive_model.rhs` |
| B-04 | `A1(y) = Alf_d*MSchkv*g*sign(y6)` | 4066 | `archive_model.rhs` |
| B-05 | `A2(y) = Alf_ch*MSkv*g*sign(y7)` | 4067 | `archive_model.rhs` |
| B-06 | `A3(y) = Alf_k*Mv*g*sign(y8)` | 4068 | `archive_model.rhs` |
| B-07 | `B1 = MprKCH/Rk − A1` | 4081 | `archive_model.rhs` |
| B-08 | `B2 = FprCHK − A2` | 4082 | `archive_model.rhs` |
| B-09 | `B3 = FprKV − A3` | 4083 | `archive_model.rhs` |
| B-10 | `BZ2 = (B1 − B2)/Mch` | 4084 | `archive_model.rhs` |
| B-11 | `BZ3 = −B1/Mch + Mchk*B2 − B3/Mkuz` | 4085 | `archive_model.rhs` |
| B-12 | `BZ4 = −B2/Mkuz + Mkv*B3` | 4086 | `archive_model.rhs` |
| B-13 | `D(t,y) = (y5,y6,y7,y8,(Mtag − MprKCH)/Itr, BZ2, BZ3, BZ4)` | 1725 | `archive_model.rhs` |
| B-14 | `D2` — braking, `Mgalm` replaces `Mtag` | 4406 | `rhs(..., braking=True)` |
| B-15 | `Mgal(t) = atan((t − T1_kin)*10)/1.55*Mgal_max` | 3804 | `ArchiveParams.Mgal` |
| B-16 | `Mgalm(t,y) = Mgal(t)*sign(y5)` | 3805 | `archive_model.rhs` |
| B-17 | `Mtag(y) = interp(cspline(om1,mt1), om1, mt1, y5)` | 2106 | `ArchiveParams.Mtag` |
| B-18 | `om1 = v1/Rk`, `mt1 = pt1*Rk` | 2103–2104 | `archive_model.V1, PT1` |

**Equivalence of B-10…B-12 to the body equations.** Substituting
`Mchk = 1/Mch + 1/Mkuz` and `Mkv = 1/Mkuz + 1/Mv` shows that BZ2–BZ4 are exactly the
relative accelerations of

```
Mch  * a_ch                      = B1 - B2
Mkuz * (a_ch + ddu1)             = B2 - B3
Mv   * (a_ch + ddu1 + ddu2)      = B3
```

Summing leaves `B1 = MprKCH/Rk − A1`, i.e. the internal friction pairs A2 and A3 cancel.
Verified numerically by `tests/test_physics.py::test_internal_forces_cancel`.

---

## C. What is **not** in the archive

| Item | Status |
|---|---|
| viscous damping in any coupling | absent; couplings are purely elastic |
| tangential wheel damping `D_kch` | absent from `MprKCH` |
| tractor as a separate mass | absent; the vehicle is a single three-body chain |
| first-gear traction characteristic | absent; see below |

The `D_kch` term used in an earlier version of this work is therefore **removed**. With
the archived parameters the model reproduces the published integral results to 1.14 %
without it, so no effective dissipation term is required and none is used.

---

## D. Cross-check against the vehicle data sheet

`ПАРАМЕТРИ ZiL-131.doc` gives kerb mass 6460 kg, driving-wheel radius 0.530 m and
reduced driveline inertia 1767 kg·m².

* wheel radius: **0.530 m in both** — independent confirmation of A-01;
* mass: `Mch + Mkuz = 6600 kg` against a kerb mass of 6460 kg — consistent to 2.2 %;
* inertia: 1767 kg·m² against `Itr = 2025 kg·m²` in the worksheet.

The worksheet value is retained, because it is the one that produced the published
results; the data-sheet value is recorded here as a plausibility check only. Mixing the
two would break reproduction.

---

## E. Other worksheets — do not mix

| Worksheet | `Mch` | `Mkuz` | cargo | `Itr` | Use |
|---|---|---|---|---|---|
| ETAL | 5900 | 700 | `Mv` 5000 | 2025 | **authoritative** |
| `Ch.P.5.2.2` | 5565 | 2700 | `Mvant` 1500 | 3746 | different case; not used |
| `Друга тест задача` | – | – | – | `IK` 1780 | simplified test case; not used |

These are distinct computational cases, not revisions of one model. Taking a convenient
value from one and another from a second would produce a model that reproduces neither.

---

## F. Unresolved

**The published Table 1 peak forces are not reproducible.** The source tabulates peak
longitudinal forces for start-off *in first gear* at cargo masses 1250–5000 kg. The
archive-locked model exceeds them by 48–58 %. The archive contains no first-gear
traction characteristic, so the scenario cannot be reconstructed. This is reported in
the manuscript as a limitation; those values are not used in any comparison and no
parameter was adjusted to approach them.
