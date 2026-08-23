# Experimental data — NOT AVAILABLE

No physical measurements have been recorded. This directory is empty by design.
The validation stage reports `NOT_AVAILABLE`; it does **not** substitute synthetic data.

## Expected CSV schema (one file per run)

File name: `run_<run_id>_rep_<replicate>.csv`

| Column | Unit | Required | Notes |
|---|---|---|---|
| `t` | s | yes | monotone, from the raw logger clock |
| `v_fifth_wheel` | m/s | yes | fifth-wheel speedometer |
| `odometer` | m | yes | fifth-wheel odometer |
| `a_x_frame` | m/s^2 | yes | frame longitudinal acceleration |
| `a_y_frame` | m/s^2 | yes | frame lateral acceleration |
| `a_z_frame` | m/s^2 | yes | frame vertical acceleration |
| `a_x_accel2` | m/s^2 | no | second accelerometer, longitudinal |
| `u1_long` | m | yes | front tank, longitudinal, relative to frame |
| `u1_vert` | m | yes | front tank, vertical, relative to frame |
| `u2_long` | m | yes | rear tank, longitudinal, relative to frame |
| `u2_vert` | m | yes | rear tank, vertical, relative to frame |
| `engine_rpm` | rpm | no | if available from the CAN bus |
| `engine_torque` | N m | no | if available |
| `engine_power` | W | no | if available |

Accompanying metadata file `run_<run_id>_rep_<replicate>.yaml`:

```yaml
run_id: 14
replicate: 2
randomized_position: 7
front_fill_pct: 50
rear_fill_pct: 50
front_tank_mass_kg:      # measured, not nominal
rear_tank_mass_kg:       # measured, not nominal
speed_factor_mps: 4.5
ambient_temperature_C:
surface_condition:
tyre_pressure_kPa:
operator:
timestamp:
notes:
```

## Processing rules

* Raw data are never overwritten; all corrections are applied downstream.
* Zero-offset correction uses the quiescent window before the manoeuvre.
* Low-pass filtering is configurable and is reported with every result.
* Missing values, outliers and the estimated sampling rate are reported, not silently
  repaired.
* For each channel report peak, RMS and the 95th percentile; for displacement also the
  peak positive, peak negative and peak absolute values.
* Across replicates report mean, standard deviation, coefficient of variation and the
  95 % confidence interval.

See `docs/EXPERIMENT_PROTOCOL.md` for the run protocol.
