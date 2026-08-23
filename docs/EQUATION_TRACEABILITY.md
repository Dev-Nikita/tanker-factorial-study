# Equation traceability

Every implemented equation, its physical meaning, its source, what was changed during
implementation, and the units. Source codes as in `MISSING_INPUTS.md`
(**SRC-A** = original tanker manuscript, **SRC-B** = 3³ factorial brief,
**SRC-C** = `red_hnups,+13.pdf`, **SRC-D** = *Теорія автомобіля* lecture notes).

---

## E-01 Tyre twist and longitudinal wheel force

**Implemented** (`model.py::rhs_rigid`, `rhs_compliant`)

    θ  = φ_k − X_ch / r_k                       [rad]
    F_kch = (C_kch·θ + D_kch·θ̇) / r_k           [N]

**Meaning.** The driving wheels have tangential (circumferential) compliance; the
twist between the wheel angle and the travelled distance generates the longitudinal
force applied to the chassis.

**Source.** SRC-A eq. (2), first and second lines — recovered from the rendered
equation object (`word/media/image39.wmf`), which reads
`I_k d²φ_k/dt² = M_k − C_kch(φ_k − X_ch/r_k)` and
`m_c d²X_ch/dt² = C_kch(φ_k − X_ch/r_k)(1/r_k) − f_d·m_c·g·sign(dX_ch/dt)`.
Consistent with assumption 3 of SRC-A ("wheels have tangential compliance depending
on the twist angle").

**Change.** The damping term `D_kch·θ̇` is **added**; SRC-A contains only `C_kch`.
Rationale and consequences: `MISSING_INPUTS.md` M-02/P-4.

**Units.** `C_kch` [N·m/rad], `D_kch` [N·m·s/rad], `r_k` [m], `F_kch` [N].

---

## E-02 Rigid baseline

**Implemented** (`model.py::rhs_rigid`)

    I_k·φ̈_k  = M_k(φ̇_k) − M_gal(t)·sgn(φ̇_k) − (C_kch·θ + D_kch·θ̇)
    m_c·Ẍ_ch = F_kch − f_d·m_c·g·sgn(Ẋ_ch)

**Meaning.** SRC-A scheme "tractor engine – tanker chassis together with the cargo":
two generalised coordinates, the whole vehicle mass moving as one body.

**Source.** SRC-A eq. (2), transcribed literally except for E-01.

**Change.** `M_gal(t)` (brake torque) is applied in the braking scenario; SRC-A states
that in eqs. (1)–(2) the engine torque `M_k` is *replaced* by the brake torque during
stage 2, which is what the `braking=True` branch does. `sgn` is regularised as
`tanh(v/ε)`, ε = 10⁻³ m/s.

**Units.** `I_k` [kg·m²], `M_k`, `M_gal` [N·m], `m_c` [kg], `f_d` [-], `g` [m/s²].

---

## E-03 Compliant–damped coupling forces

**Implemented** (`model.py::rhs_compliant`)

    F₁ = C_chv1·u₁ + α_v1·m₁·u̇₁      R₁ = f_ch·m₁·g·sgn(u̇₁)
    F₂ = C_v1v2·u₂ + α_v2·m₂·u̇₂      R₂ = f_k ·m₂·g·sgn(u̇₂)

**Meaning.** Spring plus parallel damper between chassis and cargo unit 1, and between
units 1 and 2; Coulomb resistance of the roller guides proportional to the normal load.

**Source.** SRC-A Fig. 2 (`word/media/image26.jpeg`) defines `C_chv1`, `C_v1v2`,
`α_v1`, `α_v2` and the relative coordinates `X_v1`, `X_v2` in the chassis-fixed frame
O₁X₁Z₁. SRC-A text: "α — resistance coefficients **referred to unit cargo mass**,
depending on the velocity of the first cargo and on the relative velocity of the
cargoes respectively". Hence the damping force is `α·m·u̇`, not `α·u̇`.
Guide-resistance coefficients `f_ch`, `f_k` and their range 0.02…0.40 are from SRC-A.

**Change.** Spring free lengths `l_chV1`, `l_V1V2` are absorbed into the relative
coordinates (`u = X_rel − l`), so that `u = 0` is the unstretched state. This is a
change of origin only and does not affect the dynamics.

**Units.** `C` [N/m], `α` [N·s·m⁻¹·kg⁻¹], `u` [m], `F` [N].

---

## E-04 Compliant configuration, equations of motion

**Implemented** (`model.py::rhs_compliant`, solved for the accelerations each step)

    m₀·Ẍ_ch                 = F_kch − R_road + F₁
    m₁·(Ẍ_ch + ü₁)          = −F₁ − R₁ + F₂
    m₂·(Ẍ_ch + ü₁ + ü₂)     = −F₂ − R₂
    R_road = f_d·(m₀ + m₁ + m₂)·g·sgn(Ẋ_ch),   m₀ = m_tr + m_ch

**Meaning.** SRC-A scheme "engine – chassis – series-connected tanks (cargo)"
(Fig. 2). Absolute accelerations of the units are obtained from the chassis
acceleration plus the relative ones, which is what makes the 3×3 linear solve
necessary.

**Source.** SRC-A eq. (1) *in intent*. See the caveat below.

**⚠ Caveat (also MISSING_INPUTS M-04) — REVISED after author clarification.**

An earlier version of this file stated that the display object carrying eq. (1) in the
SRC-A `.docx` (`word/media/image35.wmf`) renders as a *spatial* model "pasted from a
dissertation". **That statement was wrong and has been withdrawn.** The author of the
source model has confirmed that the original formulation is a **planar model in the
vertical longitudinal plane**: the coordinates `X`, `Z` and the rotation `φ` are the
longitudinal displacement, the vertical displacement and the pitch angle of a planar
body. There is no lateral axis `y`, so the presence of `Z` and `φ` does not indicate a
three-dimensional model.

What remains true, and is the actual relationship between E-04 and the source, is this:
the source formulation retains vertical and pitch degrees of freedom, whereas E-04 keeps
only the longitudinal ones. E-04 is therefore

    a REDUCED-ORDER LONGITUDINAL MODEL derived from the original planar
    (longitudinal-vertical-pitch) vehicle formulation,

not a literal transcription of eq. (1). The reduction is justified for the manoeuvres
studied here — straight-line start-off and braking on a horizontal road, where the
vertical and pitch responses do not couple into the longitudinal traction path at first
order — and it is declared as assumption 4 of the source itself ("the computational
scheme is treated as planar"). The manuscript must describe it in exactly these terms.

The original Mathcad worksheets are still required to confirm the transcription
term by term.

**Verification that the transcription is consistent with eq. (2).** Test
`test_stiff_limit_approaches_rigid` sets `C_chv1 = C_v1v2 = 5×10⁹ N/m` and the guide
friction to zero; the compliant model then reproduces the rigid-baseline peak wheel
force to within 2 %, as it must.

---

## E-05 Cargo mass model

**Implemented** (`config.py::unit_masses`)

    m_i = m_k + m_tank,dry + ρ·V·λ_i,   λ_i ∈ {0, 0.5, 1}

**Source.** SRC-B specifies the fill levels 0/50/100 % and the tank geometry
(50 L, ⌀370 mm, 500 mm). SRC-A gives the trolley mass `m_k = 700` kg and the cargo
range 1250…5000 kg.

**Change.** None. For the full-scale configuration `m_tank,dry` is absorbed into
`m_k` because SRC-A does not separate them.

**Carriage (bogie).** SRC-A Table 1 lists three forces (`F_kch`, `F_chk`, `F_chv`),
implying a chassis -> carriage -> cargo chain, while Fig. 2 shows chassis -> tank 1 ->
tank 2. The author has confirmed that **the carriage is omitted from Fig. 2 for clarity
only**. Each movable unit therefore lumps carriage, tank and liquid:

    m_vi = m_carriage,i + m_tank,i + m_liquid,i

which is what `unit_masses` computes. The earlier note describing this as an unresolved
ambiguity has been withdrawn.

**Units.** [kg], ρ [kg/m³], V [m³].

---

## E-06 Powertrain characteristic

**Implemented** (`model.py::Powertrain.torque`)

    ω ≤ ω_p          : M = M₀ + (M_peak − M₀)(3s² − 2s³),  s = ω/ω_p
    ω_p < ω ≤ ω_e    : M = N_d / ω,        N_d = M_peak·ω_p
    ω_e < ω < ω_max  : M = (N_d/ω_e)·(ω_max − ω)/(ω_max − ω_e)
    ω ≥ ω_max        : M = 0

with `M₀ = 1.50×10⁴` N·m, `ω_p = 2.6` rad/s, `M_peak = 3.55×10⁴` N·m,
`ω_e = 22` rad/s, `ω_max = 37` rad/s.

**Source.** SRC-A Fig. 4 (`word/media/image41.jpeg`), digitised, and SRC-A eq. (3):
`N_d = M_k · dφ/dt ≅ const`. SRC-A text: "over the angular-velocity range
2.6–22 rad/s the engine characteristic is almost ideal". SRC-A assumption 6 states
the same.

**Change.** The rising branch below `ω_p` is a smoothstep interpolation between the
two digitised end points; SRC-A gives only the curve. The cut-off branch above
`ω_e` is linear to the digitised zero-crossing.

**Units.** [N·m], [rad/s], `N_d` [W].

---

## E-07 Brake torque law

**Implemented** (`model.py::Brakes.torque`)

    M_gal(t) = s · [ M_∞ − ΔM · exp(−(t − t_b)/τ) ],  t ≥ t_b

with `M_∞ = 1750` N·m, `ΔM = 1200` N·m, `τ = 0.8` s, `t_b = 600` s.

**Source.** SRC-A Fig. 5 (`word/media/image48.jpeg`), digitised (≈550 N·m at
t = 600 s rising to ≈1700 N·m at t = 610 s).

**Change.** The scale factor `s = 3.9652` is **added** and identified from the
published braking time of 25.054 s; it is interpreted as the equivalent number of
braked wheel sets, since Fig. 5 evidently gives the torque of one set. With `s = 1`
the simulated stopping time is 49 s. With the identified `s` the stopping distance is
166.3 m against the published 163.6 m (+1.7 %), which was **not** a fitting target and
therefore serves as an independent check.

---

## E-08 Response-surface model

**Implemented** (`regression.py`)

    Y = β₀ + Σβ_i x_i + Σβ_ij x_i x_j + Σβ_ii x_i² + ε      (10 coefficients)

**Source.** Standard second-order response-surface methodology; SRC-B asks for "a
polynomial of the full factorial experiment using a three-level full-factorial plan
(3³ = 27)" for a qualitative assessment of nonlinearity.

**Change.** SRC-B could be read as requesting the full 27-term tensor polynomial.
That model is saturated (27 coefficients, 27 design points, zero residual degrees of
freedom) and supports no inference. It is fitted in
`04_fit_response_surfaces.py::saturated_check` and reported only as an interpolant,
labelled as such. The 10-coefficient model is used for every statistical statement.

---

## E-09 Energy functional used in the tests

**Implemented** (`tests/test_physics.py::_mech_energy`)

    E = ½(m₀Ẋ_ch² + m₁v₁² + m₂v₂² + I_k φ̇²) + ½(C_chv1 u₁² + C_v1v2 u₂²)
    v₁ = Ẋ_ch + u̇₁,  v₂ = Ẋ_ch + u̇₁ + u̇₂

Used to verify conservation without dissipation (drift < 10⁻⁵ relative) and monotone
decay with damping, i.e. the Lyapunov argument stated in the manuscript.

---

## Symbols used but not appearing in SRC-A

| Symbol | Meaning | Status |
|---|---|---|
| `D_kch` | tangential wheel damping | added term, E-01 |
| `m_tr` | tractor mass | identified, P-1 |
| `r_k` | dynamic wheel radius | identified, P-2 |
| `I_k` | driveline inertia at wheel axes | identified, P-3 |
| `s` | braked wheel sets | identified, E-07 |
| `u₁, u₂` | relative coordinates with free length removed | change of origin, E-03 |
| `ε` | Coulomb regularisation velocity | numerical, 10⁻³ m/s |
