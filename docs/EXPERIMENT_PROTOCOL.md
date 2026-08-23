# Physical test protocol (planned campaign)

**Status: NOT EXECUTED.** No sensor data exist. This document specifies the campaign so
that it can be run against a pre-registered design, and defines the data format the
analysis pipeline expects.

---

## 1. Rig

From SRC-B, verbatim in substance:

- two tanks, nominal volume 50 L each, ⌀370 mm, length 500 mm;
- tanks mounted on longitudinal guides on the semitrailer frame, free to travel along it;
- the front tank connected to the frame by a spring;
- the two tanks connected to each other by a spring;
- displacement transducers between frame and each tank, measuring **longitudinal and
  vertical** displacement;
- triaxial vibration accelerometers on the frame (longitudinal, lateral, vertical);
- an additional accelerometer on the frame;
- a fifth wheel with odometer and speedometer.

**To be measured and recorded before the campaign** (currently placeholders, see
`MISSING_INPUTS.md` M-05): dry tank mass, roller-carriage mass, spring rates of both
springs, damping (if dampers are fitted), guide friction, allowable travel, liquid type
and density, and the tractor data (kerb mass, tyre size, first-gear ratio).

---

## 2. Design

Three factors, three levels, full factorial, 27 unique combinations:

| Factor | Quantity | −1 | 0 | +1 |
|---|---|---|---|---|
| x₁ | front-tank fill | 0 % | 50 % | 100 % |
| x₂ | rear-tank fill | 0 % | 50 % | 100 % |
| x₃ | speed factor | 4.0 m/s | 4.5 m/s | 5.0 m/s |

The design matrix, including the randomised execution order generated from seed 42, is
`results/tables/factorial_design.csv` (columns `run_id`, `randomized_run_id`).

**⚠ Before the first run, fix the meaning of x₃** (`MISSING_INPUTS.md` M-01). The
simulation assumes *target speed of the start-off manoeuvre*. If the intended meaning
is an initial speed, the manoeuvre is a different one and the simulation must be rerun
before any comparison is made.

---

## 3. Replication

**Preferred: 27 × 3 = 81 runs**, three replicates of every combination.

Replication is what supplies the pure-error estimate. Without it there is no way to
separate lack of fit of the response surface from measurement scatter, and no honest
confidence interval on any measured effect. A deterministic simulation cannot supply
this — identical inputs return identical outputs — which is why the computational study
uses Monte-Carlo propagation instead and why the two must never be conflated.

**Minimum acceptable: 27 + 5 = 32 runs** — all 27 unique combinations plus five
replicated centre points at 50 % / 50 % / 4.5 m·s⁻¹. Five centre points give a
usable pure-error estimate with 4 degrees of freedom and detect curvature, but give no
lack-of-fit information away from the centre.

---

## 4. Randomisation

Execute in `randomized_run_id` order, not in `run_id` order.

Over a session, tyre temperature, brake temperature, surface condition, ambient
temperature and battery state all drift. Running the design in standard order aliases
that drift with the factor effects — most damagingly with x₃, which would otherwise be
executed as three consecutive blocks. Randomisation converts a systematic bias into
noise that replication can then quantify.

Re-randomise between replicate sets (three independent permutations, seeds 42, 43, 44).

---

## 5. Per-run procedure

1. Set fill levels by mass, not by volume; record the actual mass of each tank
   (±0.5 kg) rather than the nominal one.
2. Return both carriages to their reference position; record the reference reading of
   both displacement transducers.
3. Bring the vehicle to the defined initial state and let the transducers settle
   (≥ 5 s of quiescent signal must be recorded before the manoeuvre).
4. Execute the manoeuvre; log continuously from at least 2 s before to at least 10 s
   after it.
5. Record ambient temperature, surface condition and tyre pressures in the run log.
6. Do not adjust filtering, gain or zero between runs of the same replicate set.

Sampling rate: at least 500 Hz on the accelerometers (the tyre-twist transient of the
simulated system sits near 5 Hz, and jerk requires headroom for differentiation);
at least 200 Hz on the displacement channels.

---

## 6. Data format

Deposit one CSV per run in `data/experimental/` following the schema in
`data/experimental/README.md`. Raw data must be preserved unmodified; all filtering is
applied downstream and is configurable.

---

## 7. Comparison with the simulation

Two cautions.

**Scale.** The rig tanks hold 50 kg of water against a multi-tonne frame, whereas the
full-scale model carries 5000 kg per tank. The mechanism studied here depends on the
separation between the coupling frequency √(C/m) and the tyre-twist transient, so it
does not survive scaling unchanged. Running the same design at rig scale in the
repository returns a mean peak-force reduction of 0.86 % against 9.18 % at full scale.
A similitude analysis is required before the rig is used to validate the full-scale
claim; the rig validates the **model**, not the magnitude of the benefit.

**Metrics.** Compare per channel using RMSE, NRMSE, MAE, R² over the transient window,
plus the peak error

    E_peak = |F_exp,max − F_sim,max| / F_exp,max × 100 %

Report these per run and aggregated. Do not tune model parameters to the rig data and
then present the agreement as validation; if parameters are re-identified, say so and
hold out runs for testing.
