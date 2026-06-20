# Ireland Dino Finder

This started as a joke about finding dinosaur bones with sound waves in Ireland. The useful version is an R&D sandbox for weird wave math: simulate sound or pressure waves through 2D and eventually 3D material slices, with special interest in fractional differential equations whose fractional order varies over space by material.

The current Rust code is a baseline toy with useful scaffolding: a 2D scalar acoustic wave equation with Ricker source, sponge boundary damping, per-cell wave speed, CSV/HDF5 output, and fast diagnostic tests. The project goal is to grow that into a simulation and visualization pipeline for heterogeneous, dispersive, lossy media.

## Current State

- `src/main.rs`: thin demo entry point.
- `src/grid.rs`, `src/materials.rs`, `src/source.rs`, `src/model.rs`, `src/solver.rs`, `src/output.rs`: simulation modules.
- `src/hdf5_helper.rs`: HDF5 time-series writer used by the active solver path.
- `python/csv_vis.py`: quick matplotlib viewer for `output/pressure_final.csv`.
- `python/h5py_reader.py`: helper for reading future HDF5 wave cubes.
- `python/h5_to_mp4.py`: render `output/wavefield.h5` to `output/wavefield.mp4`.
- `data/material_properties.csv`: seed material properties with citation keys, source URLs, and uncertainty/range columns.
- `data/material_map_example.csv`: tiny material-id map showing the CSV mask format.
- `source/pdf/`: local source-paper cache. PDFs are ignored by git.

Run the current demo:

```sh
cargo run
python python/csv_vis.py
python python/h5_to_mp4.py
```

Run visual sanity scenarios:

```sh
cargo run -- --scenario visual-all
ffmpeg -y -framerate 30 -i output/visual/visual_lossless_uniform_frames/frame_%06d.ppm -pix_fmt yuv420p output/visual/visual_lossless_uniform.mp4
ffmpeg -y -framerate 30 -i output/visual/visual_lossless_slow_block_frames/frame_%06d.ppm -pix_fmt yuv420p output/visual/visual_lossless_slow_block.mp4
ffmpeg -y -framerate 30 -i output/visual/visual_biot_sand_target_frames/frame_%06d.ppm -pix_fmt yuv420p output/visual/visual_biot_sand_target.mp4
```

If your Python environment has `h5py`, you can also render directly from HDF5:

```sh
python python/h5_to_mp4.py output/visual/visual_lossless_slow_block.h5 output/visual/visual_lossless_slow_block.mp4
```

Run the fast baseline checks:

```sh
cargo test
```

## Research Goal

Model wave propagation through a spatially heterogeneous slice:

- water, sand, sediment, rock, bone-like inclusions, and engineered phantom materials;
- space-varying wave speed, density, attenuation, and fractional order;
- video output of wavefields, residuals, energy loss, and model comparisons;
- baseline models side-by-side with fractional models;
- eventual bench experiment: a fish tank filled with sand or glass beads, water saturation controlled as well as practical, and a bone or bone-like target embedded in the middle.

The Substack angle is the trail of failed and partially successful models: what each model predicts, what frequency-dependent behavior it cannot explain, and whether fractional/variable-order models are a useful compact language or just an attractive nuisance.

## Working Hypothesis

Many real materials do not behave like the clean acoustic wave equation. Measured attenuation and phase velocity often depend on frequency. In tissues, sediments, and porous materials this dependence is commonly modeled as a power law:

```text
attenuation(f) = a0 * f^y
```

where `y` is often non-integer and material dependent. A single integer-order damping term, a single relaxation time, or a lossless wave equation cannot generally match that behavior over a wide band. The question for this project is whether a spatially varying fractional order, for example `alpha = alpha(x, z)`, gives a useful simulation knob for material-dependent loss and dispersion.

## Candidate Fractional Model

A fractional Zener-style constitutive relation is the main target because it is connected to relaxation physics rather than being only curve-fitting:

```text
sigma + tau_sigma^alpha D_t^alpha sigma
  = E0 * (epsilon + tau_epsilon^beta D_t^beta epsilon)
```

For acoustic pressure-like simulations, one commonly sees wave-equation forms such as:

```text
laplacian(u)
  - (1 / c0^2) d_t^2 u
  + tau_sigma^alpha D_t^alpha laplacian(u)
  - (tau_epsilon^beta / c0^2) D_t^(beta + 2) u
  = s
```

The R&D extension is to let parameters vary spatially:

```text
c0 = c0(x, z)
rho = rho(x, z)
alpha = alpha(x, z)
beta = beta(x, z)
tau_sigma = tau_sigma(x, z)
tau_epsilon = tau_epsilon(x, z)
```

That is mathematically and numerically awkward. The first practical implementation should probably be a 2D material map with piecewise-constant parameters, where each material owns a fixed fractional order and relaxation parameters.

Not every exponent is physically admissible. The Caputo power-law memory kernel is completely monotone, hence a positive superposition of decaying relaxation modes, but a full attenuation/dispersion model has extra passivity constraints. In Hanyga's power-law attenuation example `beta(p) = C p^alpha`, complete monotonicity of the associated relaxation modulus requires roughly `1/2 <= alpha <= 1`. Exponents outside the relevant passive range can still be fit, but they should be flagged per material rather than treated as automatically physical.

## Baseline Models

### 1. Lossless Scalar Acoustic Wave

```text
d_t^2 p = c(x)^2 * laplacian(p) + s
```

Variable-density form:

```text
(1 / K(x)) d_t^2 p = div((1 / rho(x)) grad(p)) + s
```

Strengths: simple, fast, easy to validate; good for first arrival times, basic reflection/refraction, and code correctness.

Limits: no intrinsic attenuation, no frequency-dependent phase velocity, cannot reproduce power-law attenuation; interfaces look too clean unless scattering is explicitly resolved.

Implementation status: implemented as `WaveModel::LosslessAcoustic` in `src/model.rs`. The active solver supports homogeneous and heterogeneous per-cell wave speed maps.

### 2. Simple Damped Acoustic Wave

```text
d_t^2 p + gamma(x) d_t p = c(x)^2 * laplacian(p) + s
```

Strengths: easy first lossy model; useful for checking that energy decays and boundaries behave.

Limits: damping is not a good broadband material law; one coefficient cannot match both amplitude loss and dispersion; measured attenuation is usually closer to `f^y` than to a viscous term.

Implementation status: implemented as `WaveModel::LinearDampedAcoustic { gamma }` in `src/model.rs`.

### 3. Standard Linear Solid / Zener Relaxation

One relaxation mechanism can be represented by a frequency-dependent modulus:

```text
M(omega) = M_inf * (1 + i * omega * tau_epsilon)
                 / (1 + i * omega * tau_sigma)
```

Strengths: physically interpretable spring-dashpot baseline; causal dispersion and attenuation; extendable to multiple relaxation mechanisms.

Limits: a single relaxation time fits a limited band; many mechanisms fit broader data but become parameter-heavy; does not by itself explain non-integer broadband power laws compactly.

Implementation status: reduced baseline implemented as `WaveModel::StandardLinearSolid { damping_gamma, relaxation_time_s, relaxation_strength }` in `src/model.rs`. This version uses a one-memory relaxation of the Laplacian/strain term plus damping, so it is closer to a Zener-style relaxation baseline than the first velocity-memory proxy, but it is still not a fully calibrated constitutive solver.

### 4. Fractional Constant-Q / Kjartansson-Style Viscoacoustic Model

A common seismic approximation treats quality factor `Q` as roughly frequency independent over a band:

```text
Q^-1 = energy_lost_per_cycle / (2 * pi * stored_energy)
attenuation(omega) approximately proportional to omega / (2 * Q * c)
```

Strengths: widely used in seismic modeling; better than lossless acoustics for broad attenuation; good comparison point for inverse-Q processing; the Kjartansson idealization is a fractional-law baseline rather than a plain integer-order damper.

Limits: constant `Q` implies a constrained frequency law; low- and high-frequency behavior need care for causality; not enough when measured exponents differ strongly by material.

Implementation status: reduced fractional baseline implemented as `WaveModel::FractionalConstantQ { q, reference_freq_hz, dispersion_strength }` in `src/model.rs`. This version maps `Q` at a reference frequency to damping and adds a causal one-relaxation dispersion proxy around that frequency, so it is still band-limited rather than a full Kjartansson implementation.

### 5. Biot Poroelastic Model

For saturated porous media, Biot theory models coupled solid-frame and pore-fluid motion. A schematic form is:

```text
solid momentum: rho_11 d_t^2 u + rho_12 d_t^2 U = div(sigma)
fluid momentum: rho_12 d_t^2 u + rho_22 d_t^2 U + b d_t(U - u) = -grad(p_f)
```

where `u` is solid displacement, `U` is fluid displacement, `p_f` is pore pressure, and `b` depends on permeability and fluid viscosity.

Strengths: physics-informed for water-saturated sand and sediments; predicts multiple compressional modes and strong attenuation; directly relevant to the fish-tank experiment.

Limits: many hard-to-measure parameters (permeability, tortuosity, frame moduli, pore geometry); much heavier to implement than scalar acoustics; published measurements in water-saturated granular media still show frequency dependence that simpler Biot-derived models may only match qualitatively.

Implementation status: reduced baseline implemented as `WaveModel::ReducedBiotPoroelastic { drag_gamma, relaxation_time_s, pore_coupling }` in `src/model.rs`. This first version is a pore-drag memory proxy for smoke testing, not full coupled solid/fluid Biot elastodynamics.

Current visual scenarios:

- `visual_lossless_uniform`: basic circular wavefront sanity check.
- `visual_lossless_slow_block`: wavefront delay/refraction through a slow rectangular block.
- `visual_lossless_fast_target`: sand-like slow region with a fast circular target.
- `visual_damped_slow_block`: same slow block with simple damping.
- `visual_sls_slow_block`: slow block with the reduced SLS/Zener baseline.
- `visual_fractional_constant_q_slow_block`: slow block with the reduced fractional constant-Q baseline.
- `visual_biot_sand_target`: slow sand-like region with fast target and reduced Biot/EDFM-style drag.

Each scenario writes:

- `output/visual/<name>.h5`: HDF5 wavefield frames.
- `output/visual/<name>.csv`: final pressure field.
- `output/visual/<name>_frames/frame_*.ppm`: dependency-free image frames for `ffmpeg`.

## Why Frequency Dependence Matters

The baseline models should not only be compared by "does the wave look nice." They should be compared against band-dependent measurements:

- arrival time or phase velocity as a function of frequency;
- amplitude attenuation as a function of frequency;
- reflection/transmission coefficients at material interfaces;
- recovered target contrast under matched processing.

Attenuation and dispersion are not two independent knobs. Causality forces them into a Hilbert-transform pair (the Kramers-Kronig relations): fix the attenuation law and the dispersion is essentially determined. A model that matches measured attenuation but predicts the wrong phase velocity is usually violating causality, which is why causal power-law wave equations (Szabo) are built the way they are, and why the comparisons above should check attenuation *and* dispersion together rather than tuning each separately.

Argo et al. measured water-saturated glass beads from 300 kHz to 800 kHz and reported porosity-dependent sound speed plus negative dispersion above about 550 kHz. That is exactly the kind of behavior to use as a sanity check: a model can match one frequency and still miss the trend.

## Fast Baseline Test

The first automated benchmark is intentionally tiny. It is not a research-grade validation case; it is a fast smoke test that almost any reasonable wave model should pass before we trust it on heterogeneous materials.

The test geometry is a homogeneous 2D grid:

```text
grid: 61 x 61
dx = dz = 0.05
c0 = 1.0
dt = 0.015
steps = 150
source: Ricker wavelet, f0 = 5.0, centered at t = 3 / f0
probe: 1.0 distance unit to the right of the source
boundaries: no sponge damping for the test
```

Acceptance checks:

- model 1, lossless acoustics, must produce a measurable probe pulse whose peak arrives near `source_peak_time + distance / c0`;
- lossy baseline models must keep roughly the same travel time as model 1;
- lossy baseline models must reduce peak amplitude and probe-trace energy relative to model 1.

Extra debugging checks:

- impulse sources fire only on the selected step;
- a zero-source run stays exactly at rest;
- symmetric probes in a homogeneous lossless field agree;
- a heterogeneous slow-zone path delays the probe peak;
- sponge boundaries reduce late field energy;
- unstable Courant numbers are rejected before time stepping;
- probe traces and final fields stay finite.
- the material-property CSV parser handles quoted fields and material-id masks;
- the active solver writes final CSV and HDF5 frames;
- prototype Grünwald-Letnikov fractional weights have expected first terms.

This gives future models a minimum contract:

- do not explode numerically on a small stable grid;
- preserve first-order travel time in a homogeneous medium;
- apply loss only when the model claims to apply loss;
- produce probe traces that can later be compared automatically.

The current tests live with the solver/model/material modules and run with `cargo test`.

## Data Plan

`data/material_properties.csv` is a seed file, not a truth database. Every row has:

- a material id and state;
- density and wave-speed values where available;
- derived P-wave modulus where it makes sense;
- porosity and frequency band when the value is band-specific;
- a citation key and source URL;
- notes explaining approximations.
- uncertainty/range columns for density, P-wave speed, and attenuation when a useful range is known.

Material maps can be loaded from CSV masks whose cells are material ids from `data/material_properties.csv`. See `data/material_map_example.csv` for the current simple format. Missing velocities fall back to a caller-provided speed so early experiments can still run while the dataset is incomplete.

Material maps can also be loaded from ASCII P3 PPM image masks by providing a color-to-material-id table. The loader maps each pixel color to a material id, then resolves wave speed through the same material catalog.

Immediate cleanup targets:

- replace generic cortical bone values with the actual bone, antler, or phantom used in a tank;
- add natural sand measurements, not only glass bead analogs;
- separate intrinsic material loss from scattering loss;
- add Ireland-relevant soils/rocks only when sourced from real geotechnical or geological measurements.

## Fractional-Order Interface Decision

For the first variable-order fractional implementation, treat material parameters as piecewise constant per grid cell. A material id owns `alpha`, relaxation times, and any history approximation parameters. At interfaces, update each cell using its own parameters and let the finite-difference stencil couple neighboring pressures.

This is deliberately conservative:

- it matches the current material-map representation;
- it avoids inventing an interpolation rule before we have validation data;
- it lets sharp interfaces represent real material boundaries;
- smoothing or harmonic/interface averaging can be added later if tests show numerical artifacts.

The practical consequence is that fractional history storage should eventually be grouped by material region or by quantized `alpha`, not by arbitrary floating-point values at every cell.

See `docs/fractional_numerics_pitfalls.md` for the current theory-side warning list: short-memory truncation, SOE approximations, nonsmooth-data convergence loss, variable-order issues, spatial fractional boundary traps, and scheme stability, plus a pre-merge implementation checklist. See `docs/variable_space_fractional_wave_review.md` for a longer review draft aimed at the eventual blog post.

## Output Plan

Recently completed:

- store material speed maps, material-id indices, source wavelet, probe coordinates, and probe traces alongside wavefields in HDF5.

Medium-term:

- compare lossless, damped, Zener, fractional constant-Q, Biot/EDFM, and full fractional models on the same geometry;
- compute residual plots and frequency-domain probe diagnostics;
- add richer 2D heterogeneous material maps from image masks.

Long-term:

- 3D slices;
- inverse experiments: infer target location or material map from boundary sensors;
- fish-tank validation with measured source and receiver traces.

## TODO

- [x] Refactor `main.rs` into modules: grid, materials, source, solver, output.
- [x] Add a material-map loader from CSV material-id masks.
- [x] Implement HDF5 frame writing in the active solver path.
- [x] Add video rendering from HDF5.
- [x] Add baseline lossless homogeneous acoustic solver.
- [x] Add simple damped acoustic baseline.
- [x] Add reduced Zener/SLS baseline.
- [x] Add reduced fractional constant-Q baseline.
- [x] Add reduced Biot/EDFM-style baseline.
- [x] Add baseline lossless heterogeneous acoustic solver.
- [x] Improve reduced Zener/SLS proxy from velocity memory to Laplacian/strain relaxation memory.
- [x] Add causal single-relaxation dispersion proxy to fractional constant-Q baseline.
- [ ] Replace reduced Zener/SLS proxy with calibrated constitutive model.
- [ ] Replace reduced fractional constant-Q proxy with full causal Kjartansson-style implementation.
- [ ] Replace reduced Biot proxy with coupled poroelastic or EDFM implementation.
- [x] Add named visual sanity scenarios for video inspection.
- [x] Prototype fractional time derivative weights with Grünwald-Letnikov coefficients.
- [x] Decide how to handle spatially varying fractional order at material interfaces.
- [x] Add unit tests for CFL checks, source wavelet, and material-map parsing.
- [x] Add first benchmark case with analytic travel-time behavior.
- [x] Add uncertainty columns to material data.
- [x] Add image-mask material-map loader.
- [x] Store material maps, source wavelet, and probe traces in HDF5 metadata.
- [ ] Add published-reference benchmark behavior from Argo et al. or another source paper.
- [ ] Add fractional-kernel approximation tests from `docs/fractional_numerics_pitfalls.md`.

## Local Papers

The following PDFs were downloaded into `source/pdf/` for local reading and are ignored by git:

- `nasholm_holm_2012_fractional_acoustic_wave_equation.pdf`
- `nasholm_holm_2013_fractional_zener_elastic_wave_equation.pdf`
- `holm_nasholm_2013_comparison_fractional_wave_equations.pdf`
- `baker_banjai_2020_numerical_analysis_lossy_power_law_wave_equation.pdf`
- `chintada_rau_goksel_2022_spectral_ultrasound_sos_attenuation.pdf`
- `tsiklauri_2002_poroelastic_biot_slip_velocity.pdf`
- `argo_guild_wilson_2009_sound_speed_water_saturated_glass_beads.pdf`
- `hanyga_2013_viscoelastic_completely_monotonic_relaxation.pdf`
- `lewandowska_kosztolowicz_2006_short_memory_subdiffusion.pdf`
- `jin_lazarov_zhou_2015_l1_nonsmooth_data.pdf`
- `jin_lazarov_zhou_2018_nonsmooth_data_overview.pdf`
- `li_wang_xie_2019_l1_fractional_wave_nonsmooth.pdf`
- `zheng_2021_variable_order_integral_equation.pdf`
- `quan_wu_yang_2022_fast_l2_stability_soe.pdf`
- `chaudhary_diethelm_farhadi_fuchs_2025_exponential_sum_power_law.pdf`
- `lischke_pang_gulian_2018_what_is_fractional_laplacian.pdf`
- `borthagaray_leykekhman_nochetto_2020_fractional_laplacian_boundary_singularity.pdf`

## References

- [nasholm2012] S. P. Nasholm and S. Holm, "A Fractional Acoustic Wave Equation from Multiple Relaxation Loss and Conservation Laws", arXiv:1202.4251, https://arxiv.org/abs/1202.4251.
- [nasholm2013] S. P. Nasholm and S. Holm, "On a Fractional Zener Elastic Wave Equation", arXiv:1212.4024, https://arxiv.org/abs/1212.4024.
- [holm2013] S. Holm and S. P. Nasholm, "Comparison of fractional wave equations for power law attenuation in ultrasound and elastography", arXiv:1306.6507, https://arxiv.org/abs/1306.6507.
- [baker2020] K. Baker and L. Banjai, "Numerical analysis of a wave equation for lossy media obeying a frequency power law", arXiv:2012.04520, https://arxiv.org/abs/2012.04520.
- [argo2009] T. F. Argo IV, M. D. Guild, P. S. Wilson, M. Schroter, C. Radin, and H. L. Swinney, "Sound speed in water-saturated glass beads as a function of frequency and porosity", arXiv:0906.4798, https://arxiv.org/abs/0906.4798.
- [chintada2022] B. R. Chintada, R. Rau, and O. Goksel, "Spectral Ultrasound Imaging of Speed-of-Sound and Attenuation Using an Acoustic Mirror", arXiv:2201.01435, https://arxiv.org/abs/2201.01435.
- [tsiklauri2002] D. Tsiklauri, "Phenomenological model of propagation of the elastic waves in a fluid-saturated porous solid with non-zero boundary slip velocity", arXiv:physics/0201045, https://arxiv.org/abs/physics/0201045.
- [mainardi2010] F. Mainardi, "Fractional Calculus and Waves in Linear Viscoelasticity: An Introduction to Mathematical Models", Imperial College Press, London, 2010 (2nd ed., World Scientific, 2022).
- [hanyga2013] A. Hanyga, "Wave propagation in linear viscoelastic media with completely monotonic relaxation moduli", Wave Motion 50 (2013) 909-928, arXiv:1302.0402, https://arxiv.org/abs/1302.0402.
- [szabo1994] T. L. Szabo, "Time domain wave equations for lossy media obeying a frequency power law", J. Acoust. Soc. Am. 96 (1994) 491-500.
- [waters2005] K. R. Waters, J. Mobley, and J. G. Miller, "Causality-imposed (Kramers-Kronig) relationships between attenuation and dispersion", IEEE Trans. Ultrason. Ferroelectr. Freq. Control 52 (2005) 822-833.
- [water_data_and_attenuation] Water physical-property table and ultrasound attenuation summary, https://en.wikipedia.org/wiki/Water_(data_page) and https://en.wikipedia.org/wiki/Attenuation.
- [density_logging] Matrix density values used in density logging, https://en.wikipedia.org/wiki/Density_logging.
- [p_wave_table] Representative P-wave velocity ranges for common rock types, https://en.wikipedia.org/wiki/P_wave.
- [culjat2009_attenuation] M. O. Culjat, D. Goldenberg, P. Tewari, and R. S. Singh, "A Review of Tissue Substitutes for Ultrasound Imaging", Ultrasound in Medicine & Biology, cited through the attenuation table at https://en.wikipedia.org/wiki/Attenuation.
