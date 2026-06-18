# Ireland Dino Finder

This started as a joke about finding dinosaur bones with sound waves in Ireland. The useful version is an R&D sandbox for weird wave math: simulate sound or pressure waves through 2D and eventually 3D material slices, with special interest in fractional differential equations whose fractional order varies over space by material.

The current Rust code is only a baseline toy: a constant-velocity 2D scalar acoustic wave equation with a Ricker source, sponge boundary damping, and final pressure output as CSV. The project goal is to grow that into a simulation and visualization pipeline for heterogeneous, dispersive, lossy media.

## Current State

- `src/main.rs`: 2D finite-difference pressure wave demo with constant `c0`, plus fast baseline tests for README models 1 and 2.
- `src/hdf5_helper.rs`: draft HDF5 time-series writer for future video-ready outputs.
- `python/csv_vis.py`: quick matplotlib viewer for `output/pressure_final.csv`.
- `python/h5py_reader.py`: helper for reading future HDF5 wave cubes.
- `data/material_properties.csv`: seed material properties with citation keys and source URLs.
- `source/pdf/`: local source-paper cache. PDFs are ignored by git.

Run the current demo:

```sh
cargo run
python python/csv_vis.py
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

## Baseline Models

### 1. Lossless Scalar Acoustic Wave

```text
d_t^2 p = c(x)^2 * laplacian(p) + s
```

Variable-density form:

```text
(1 / K(x)) d_t^2 p = div((1 / rho(x)) grad(p)) + s
```

Upsides:

- simple, fast, easy to validate;
- useful for first arrival times and basic reflection/refraction behavior;
- good numerical baseline for code correctness.

Downsides and limitations:

- no physical attenuation;
- no frequency-dependent phase velocity;
- cannot reproduce measured power-law attenuation;
- material interfaces can look too clean unless scattering is explicitly resolved.

Implementation status: implemented as `WaveModel::LosslessAcoustic` in `src/main.rs`.

### 2. Simple Damped Acoustic Wave

```text
d_t^2 p + gamma(x) d_t p = c(x)^2 * laplacian(p) + s
```

Upsides:

- easy first lossy model;
- useful for checking that energy decays and boundaries behave.

Downsides and limitations:

- damping is not a good broadband material law;
- one coefficient cannot usually match both amplitude loss and dispersion;
- measured attenuation is often closer to `f^y` than to a simple viscous term.

Implementation status: implemented as `WaveModel::LinearDampedAcoustic { gamma }` in `src/main.rs`.

### 3. Standard Linear Solid / Zener Relaxation

One relaxation mechanism can be represented by a frequency-dependent modulus:

```text
M(omega) = M_inf * (1 + i * omega * tau_epsilon)
                 / (1 + i * omega * tau_sigma)
```

Upsides:

- physically interpretable spring-dashpot baseline;
- causal dispersion and attenuation;
- can be extended to multiple relaxation mechanisms.

Downsides and limitations:

- a single relaxation time fits a limited frequency band;
- many relaxation mechanisms can fit data but become parameter-heavy;
- does not by itself explain non-integer broadband power laws compactly.

### 4. Constant-Q / Kjartansson-Style Viscoacoustic Model

A common seismic approximation treats quality factor `Q` as roughly frequency independent over a band:

```text
Q^-1 = energy_lost_per_cycle / (2 * pi * stored_energy)
attenuation(omega) approximately proportional to omega / (2 * Q * c)
```

Upsides:

- widely used in seismic modeling;
- better than lossless acoustics for broad attenuation;
- good comparison point for inverse-Q style processing.

Downsides and limitations:

- constant `Q` implies a constrained frequency law;
- low-frequency and high-frequency behavior need care for causality;
- not enough when measured exponents differ strongly by material.

### 5. Biot Poroelastic Model

For saturated porous media, Biot theory models coupled solid-frame and pore-fluid motion. A schematic form is:

```text
solid momentum: rho_11 d_t^2 u + rho_12 d_t^2 U = div(sigma)
fluid momentum: rho_12 d_t^2 u + rho_22 d_t^2 U + b d_t(U - u) = -grad(p_f)
```

where `u` is solid displacement, `U` is fluid displacement, `p_f` is pore pressure, and `b` depends on permeability and fluid viscosity.

Upsides:

- physics-informed for water-saturated sand and sediments;
- predicts multiple compressional modes and strong attenuation;
- directly relevant to the fish-tank experiment.

Downsides and limitations:

- many hard-to-measure parameters: permeability, tortuosity, frame moduli, pore geometry;
- implementation is much heavier than scalar acoustics;
- published measurements in water-saturated granular materials still show frequency dependencies that simpler Biot-derived models may only match qualitatively.

## Why Frequency Dependence Matters

The baseline models should not only be compared by "does the wave look nice." They should be compared against band-dependent measurements:

- arrival time or phase velocity as a function of frequency;
- amplitude attenuation as a function of frequency;
- reflection/transmission coefficients at material interfaces;
- recovered target contrast under matched processing.

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
- model 2, simple damped acoustics, must keep roughly the same travel time as model 1;
- model 2 must reduce peak amplitude and probe-trace energy relative to model 1.

This gives future models a minimum contract:

- do not explode numerically on a small stable grid;
- preserve first-order travel time in a homogeneous medium;
- apply loss only when the model claims to apply loss;
- produce probe traces that can later be compared automatically.

The current tests live in `src/main.rs` and run with `cargo test`.

## Data Plan

`data/material_properties.csv` is a seed file, not a truth database. Every row has:

- a material id and state;
- density and wave-speed values where available;
- derived P-wave modulus where it makes sense;
- porosity and frequency band when the value is band-specific;
- a citation key and source URL;
- notes explaining approximations.

Immediate cleanup targets:

- replace generic cortical bone values with the actual bone, antler, or phantom used in a tank;
- add natural sand measurements, not only glass bead analogs;
- separate intrinsic material loss from scattering loss;
- store uncertainty ranges instead of single values;
- add Ireland-relevant soils/rocks only when sourced from real geotechnical or geological measurements.

## Output Plan

Near-term:

- write HDF5 wavefield frames at interval `N`;
- add a Python script that renders frames to MP4;
- store material maps alongside wavefields;
- include source wavelet and probe traces in HDF5 metadata.

Medium-term:

- compare lossless, damped, Zener, constant-Q, Biot/EDFM, and fractional models on the same geometry;
- compute residual plots and frequency-domain probe diagnostics;
- add 2D heterogeneous material maps from CSV parameters.

Long-term:

- 3D slices;
- inverse experiments: infer target location or material map from boundary sensors;
- fish-tank validation with measured source and receiver traces.

## TODO

- [ ] Refactor `main.rs` into modules: grid, materials, source, solver, output.
- [ ] Add a material-map loader from CSV or image masks.
- [ ] Implement HDF5 frame writing in the active solver path.
- [ ] Add video rendering from HDF5.
- [x] Add baseline lossless homogeneous acoustic solver.
- [x] Add simple damped acoustic baseline.
- [ ] Add baseline lossless heterogeneous acoustic solver.
- [ ] Add Zener/SLS baseline.
- [ ] Add constant-Q baseline.
- [ ] Add Biot or effective-density-fluid-model baseline for saturated granular media.
- [ ] Prototype fractional time derivative with convolution quadrature or diffusive approximation.
- [ ] Decide how to handle spatially varying fractional order at material interfaces.
- [ ] Add unit tests for CFL checks, source wavelet, and material-map parsing.
- [ ] Add benchmark cases with analytic or published reference behavior.
- [ ] Add uncertainty columns to material data.

## Local Papers

The following PDFs were downloaded into `source/pdf/` for local reading and are ignored by git:

- `nasholm_holm_2012_fractional_acoustic_wave_equation.pdf`
- `nasholm_holm_2013_fractional_zener_elastic_wave_equation.pdf`
- `holm_nasholm_2013_comparison_fractional_wave_equations.pdf`
- `baker_banjai_2020_numerical_analysis_lossy_power_law_wave_equation.pdf`
- `chintada_rau_goksel_2022_spectral_ultrasound_sos_attenuation.pdf`
- `tsiklauri_2002_poroelastic_biot_slip_velocity.pdf`
- `argo_guild_wilson_2009_sound_speed_water_saturated_glass_beads.pdf`

## References

- [nasholm2012] S. P. Nasholm and S. Holm, "A Fractional Acoustic Wave Equation from Multiple Relaxation Loss and Conservation Laws", arXiv:1202.4251, https://arxiv.org/abs/1202.4251.
- [nasholm2013] S. P. Nasholm and S. Holm, "On a Fractional Zener Elastic Wave Equation", arXiv:1212.4024, https://arxiv.org/abs/1212.4024.
- [holm2013] S. Holm and S. P. Nasholm, "Comparison of fractional wave equations for power law attenuation in ultrasound and elastography", arXiv:1306.6507, https://arxiv.org/abs/1306.6507.
- [baker2020] K. Baker and L. Banjai, "Numerical analysis of a wave equation for lossy media obeying a frequency power law", arXiv:2012.04520, https://arxiv.org/abs/2012.04520.
- [argo2009] T. F. Argo IV, M. D. Guild, P. S. Wilson, M. Schroter, C. Radin, and H. L. Swinney, "Sound speed in water-saturated glass beads as a function of frequency and porosity", arXiv:0906.4798, https://arxiv.org/abs/0906.4798.
- [chintada2022] B. R. Chintada, R. Rau, and O. Goksel, "Spectral Ultrasound Imaging of Speed-of-Sound and Attenuation Using an Acoustic Mirror", arXiv:2201.01435, https://arxiv.org/abs/2201.01435.
- [tsiklauri2002] D. Tsiklauri, "Phenomenological model of propagation of the elastic waves in a fluid-saturated porous solid with non-zero boundary slip velocity", arXiv:physics/0201045, https://arxiv.org/abs/physics/0201045.
- [water_data_and_attenuation] Water physical-property table and ultrasound attenuation summary, https://en.wikipedia.org/wiki/Water_(data_page) and https://en.wikipedia.org/wiki/Attenuation.
- [density_logging] Matrix density values used in density logging, https://en.wikipedia.org/wiki/Density_logging.
- [p_wave_table] Representative P-wave velocity ranges for common rock types, https://en.wikipedia.org/wiki/P_wave.
- [culjat2009_attenuation] M. O. Culjat, D. Goldenberg, P. Tewari, and R. S. Singh, "A Review of Tissue Substitutes for Ultrasound Imaging", Ultrasound in Medicine & Biology, cited through the attenuation table at https://en.wikipedia.org/wiki/Attenuation.
