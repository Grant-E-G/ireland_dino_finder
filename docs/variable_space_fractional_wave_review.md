# Variable-Space Fractional Wave Models: Review Notes

These are working notes for a longer blog post. They assume comfort with analysis, operators, dynamical systems, and asymptotics. Acoustics vocabulary and the numerics of fractional PDEs are introduced as needed, since those are the parts least likely to be familiar.

The motivating question is:

> Can a wave simulation use spatially varying fractional differential equations to model how sound behaves in different materials, without throwing away the very memory/dispersion effects that made fractional equations interesting?

The short answer from the sources is: maybe, but only if we are careful. Fractional models are not magic damping terms. They are nonlocal-in-time or nonlocal-in-space operators with real numerical and modeling traps. Several common approximations quietly collapse the model back into a finite-memory relaxation system, which may be useful, but is no longer the thing we thought we were testing.

A note on naming. "Variable-space fractional" is ambiguous, and the ambiguity matters. It can mean (a) a fractional derivative *in time* whose order and parameters vary with spatial position by material, or (b) a fractional Laplacian, i.e. a fractional operator *in space*. These are different models with different difficulties. This project is, for now, almost entirely about (a). Spatial fractional operators are treated separately in §5.5.

## 1. Physics Setup: What Is Being Modeled?

For a scalar acoustic approximation, the clean baseline is pressure `p(x, t)` satisfying

```text
d_t^2 p = c(x)^2 Δp + s.
```

with `c(x)` the local wave speed, `s` a source, and `Δ` the spatial Laplacian. In the ideal case this is conservative: it carries no intrinsic loss to heat, pore-fluid motion, microstructure, or sub-grid scattering.

Real materials do lose energy, and the parts that matter here are not captured by wave speed alone. The relevant vocabulary:

- attenuation: amplitude decreases as the wave propagates;
- dispersion: phase velocity depends on frequency;
- scattering: energy moves into unresolved directions/modes;
- poroelastic loss: pore fluid and solid frame move relative to one another;
- relaxation: internal degrees of freedom lag behind stress/strain.

The experimental hook is that attenuation is often fit by a power law over a frequency band:

```text
attenuation(f) = a0 f^y,
```

where `y` is generally non-integer and material-dependent. A non-integer exponent is exactly what integer-order damping cannot produce compactly, and it is the door through which fractional equations enter.

## 2. Why Fractional Derivatives Appear

The integer-order models give a small menu of frequency behaviors. A single viscous term `γ d_t p` does not give broadband power-law attenuation and dispersion. A single Zener / standard-linear-solid relaxation gives one characteristic time scale; several relaxations give several.

A fractional derivative encodes a continuum-like distribution of relaxation times. For a Caputo derivative of order `0 < α < 1`,

```text
D_t^α u(t) = 1 / Γ(1 - α) ∫_0^t (t - τ)^(-α) u'(τ) dτ.
```

The whole content is in the kernel `(t - τ)^(-α)`: it decays slowly, so the present state depends on the entire past with power-law memory. In frequency space, fractional derivatives act as powers of `iω`:

```text
D_t^α  ↔  (iω)^α,
```

which is the compact connection to non-integer power-law attenuation and dispersion.

**Complete monotonicity is the organizing fact.** The power-law kernel `t^(-α)` is completely monotone, so by Bernstein's theorem it is a positive superposition of decaying exponentials,

```text
t^(-α) = ∫_0^∞ e^(-λ t) dμ(λ),   μ ≥ 0.
```

This single property does most of the conceptual work later:

- it is why a sum-of-exponentials approximation is the natural discretization, not a hack (§5.2);
- it is why the *signs* of the SOE weights matter — positive weights keep the approximate kernel completely monotone, which is what preserves passivity;
- it is the link between "fractional memory" and "a spectrum of relaxation modes," i.e. the bridge between the fractional model and the Zener baselines.

Hanyga's analysis sharpens this for one important model class: if the power-law attenuation/dispersion function has the form `β(p) = C p^α`, complete monotonicity of the associated relaxation modulus requires `1/2 ≤ α ≤ 1`. That is a real constraint on which exponents are physically admissible for that model, not a generic restriction on every fractional kernel.

**Attenuation and dispersion are not independent.** Causality forces them into a Hilbert-transform pair (the Kramers-Kronig relations): once the attenuation law is fixed, the dispersion is essentially determined, and vice versa. The practical consequence is that a model which matches measured attenuation but predicts the wrong phase velocity is usually violating causality somewhere. This is the rigorous reason the project cannot tune loss and dispersion as two free knobs, and it is why Szabo-style causal power-law wave equations are constructed the way they are.

For materials, the physical story is usually phrased through a fractional Zener constitutive relation:

```text
σ + τ_σ^α D_t^α σ = E0 (ε + τ_ε^β D_t^β ε).
```

Stress and strain are related through a memory law: neither purely elastic nor simply viscous, with the fractional exponents tuning how the material remembers deformation.

## 3. The Space-Varying Idea

The project idea is not only "use a fractional derivative." It is to let the constitutive parameters vary by material:

```text
α = α(x),  β = β(x),  τ_σ = τ_σ(x),  τ_ε = τ_ε(x),  c = c(x),  ρ = ρ(x).
```

Water might be close to ordinary acoustics; saturated sand might have strong poroelastic loss; a bone-like inclusion a different attenuation exponent and velocity. A grid cell's material id determines its local parameters.

The conservative first choice is piecewise-constant parameters, `α(x) = α_j` on material region `j`. This matches how we already store material maps and avoids pretending we know how to physically interpolate a fractional order across a sharp boundary. The first model should therefore look like a material-indexed operator family rather than a globally smooth coefficient field:

```text
cell i has material m(i); cell i uses α_m, τ_m, c_m, attenuation params, history approximation m.
```

That has a dynamical-systems flavor: each region carries local internal memory state, and the pressure field couples those local memory systems through the spatial stencil.

One subtlety to flag up front: `τ_σ^α` carries units of `time^α`, so under `α = α(x)` the "relaxation time" has a spatially varying fractional dimension. This is harmless within a single material but makes naive cross-material comparison and interface coupling dimensionally awkward, which is a second reason to start piecewise-constant and quantized rather than smoothly varying.

## 4. Baselines

The README documents the baseline models and their implementation status in full (`WaveModel::LosslessAcoustic`, `LinearDampedAcoustic`, `StandardLinearSolid`, `ConstantQ`, `ReducedBiotPoroelastic`). The point here is not to restate them but to record *why each one is in the comparison set*. The baselines are not busywork; they are how we keep ourselves honest, because the fractional model only earns its place by predicting something they cannot.

- **Lossless acoustic.** Tests geometry, wave speed, finite-difference stability, and visual sanity. It can give arrival times and reflection/refraction; it cannot give any intrinsic, frequency-dependent, or sub-grid loss. If this fails, nothing fractional matters.
- **Simple damped acoustic.** The minimal lossy model. Useful only to confirm that energy decays and boundaries behave. The sharp test: if the fractional model does not differ from this spectrally, the fractional machinery is dead weight.
- **Standard linear solid / Zener.** One relaxation mechanism, with a characteristic time and causal attenuation/dispersion over a limited band. This baseline matters most because many fractional *approximations* collapse into exactly "some number of Zener mechanisms" (see §5.2). When that happens it is not wrong, but it changes the claim.
- **Constant-Q / Kjartansson.** Encodes broadband attenuation compactly via a quality factor, and sits between simple damping and a fractional constitutive law. Causal constant-Q is subtle: low- and high-frequency cutoffs matter.
- **Biot / poroelastic.** For saturated sand or sediment, scalar acoustics is too thin; Biot couples solid-frame and pore-fluid motion and is directly relevant to the fish-tank geometry (sand + water + a bone-like inclusion). The cost is parameter burden — permeability, tortuosity, porosity, frame moduli, pore geometry all drive the model, and a bad Biot parameter set can be less honest than a humble scalar baseline. Argo et al.'s water-saturated glass-bead measurements are the cautionary data point: even Biot-derived effective-density-fluid models can match porosity trends while missing the frequency trend (they report negative dispersion above ~550 kHz).

## 5. What the Fractional Literature Warns Us About

The strongest theme from the numerical-analysis sources is that fractional equations are fragile under approximation. The model's advantage comes from long memory and nonlocality, and practical approximations tend to attack exactly those properties.

One framing caveat before the specifics: much of the cited numerical-analysis literature concerns *subdiffusion*, `0 < α < 1`. The wave setting mixes notation: fractional Zener material exponents usually remain in `0 < α, β <= 1`, while direct fractional-wave equations may use effective orders in `1 < γ < 2` or include higher-order terms such as `D_t^(β + 2)`. The qualitative lessons transfer, but the regularity and stability details do not transfer automatically, and the relevant references differ. Jin-Lazarov-Zhou and Lewandowska-Kosztolowicz are subdiffusion results; Li-Wang-Xie is the in-regime fractional-*wave* nonsmooth-data analysis. It is worth tagging each warning by the regime it was proved in.

### 5.1 Short Memory Can Kill the Point

A fractional derivative remembers the whole past with a power-law tail. Truncating history to a fixed recent window may leave us without a fractional model in the intended sense. Lewandowska and Kosztolowicz warn against the short-memory principle for subdiffusion; the setting is diffusion, but the modeling lesson transfers: a long-memory process is not faithfully represented by an arbitrary finite window.

For this project:

- no fixed short-memory truncation as the default;
- if used, it must be labeled and tested against full-history behavior;
- visual similarity is not enough — compare probe spectra and attenuation slopes.

### 5.2 Sum-of-Exponentials Is Useful but Dangerous

The standard computational trick approximates the power-law kernel by

```text
t^(-α) ≈ Σ_j w_j exp(-λ_j t),
```

turning the history integral into recursive internal states `z_j'(t) = -λ_j z_j(t) + input(t)`. This is exactly the Bernstein representation of §2 made finite and discrete, which is why it is the natural method — and also why it is dangerous. A finite SOE *is* a finite collection of relaxation modes. With too few exponentials or a poorly chosen fit interval, we no longer have a broadband fractional law; we have a many-relaxation Zener model. That may be a perfectly good engineering model, but it is a different claim.

Because the exact kernel is completely monotone, the honest approximation should preserve that: positive weights keep the approximate kernel completely monotone and the scheme passive. Negative weights can fit the curve while quietly breaking the physics.

For this project:

- store the approximation interval `[δ, T]`, the weights, the decays, and the kernel-fit error;
- prefer / check positivity of weights, and require stability checks when the theory needs them;
- test whether the resulting attenuation still follows the intended power law over the source band.

### 5.3 High-Order Claims Often Require Smoothness We Do Not Have

Fractional evolution equations have weak singularities near `t = 0`; even with smooth forcing the solution is less regular than an integer-order PDE solution. Jin, Lazarov, and Zhou show classical convergence claims for the L1 scheme depend on restrictive smoothness, and nonsmooth data degrade the observed rate. Our simulations deliberately include Ricker pulses, sharp material interfaces, abrupt parameter jumps, and possibly non-smooth source turn-on, so a scheme advertised as `O(τ^(2-α))` may behave closer to first order in the cases we care about.

For this project:

- convergence tests must include nonsmooth pulse/interface cases, not only smooth manufactured solutions;
- startup corrections or graded time meshes may be needed;
- any blog claim about accuracy should distinguish formal order from observed order.

### 5.4 Variable Order Is Not a Text Substitution

Replacing `α` by `α(x)` looks harmless but can break proof structures that hold for constant order. Zheng's variable-order work notes that discretization coefficients can lose monotonicity, which blocks standard numerical-analysis techniques. This is directly relevant to material maps: at a boundary the memory law changes discontinuously between neighboring cells.

For this project:

- begin with piecewise-constant material orders;
- treat monotonicity/positivity of discrete weights as a diagnostic;
- add a two-material interface benchmark before any elaborate visuals;
- avoid smoothing `α(x)` unless it has physical or numerical justification.

### 5.5 Spatial Fractional Operators Are Boundary-Condition Minefields

"Space fractional" can also mean a fractional Laplacian `(-Δ)^s`. On all of `R^n` this is already nonlocal; on a bounded box, "the fractional Laplacian" is underspecified — integral, spectral, directional, and other definitions genuinely differ, and boundary behavior is part of the operator. Lischke et al. survey the definitions; Borthagaray, Leykekhman, and Nochetto show boundary-singularity effects for the integral version. For wave simulation with absorbing boundaries this is serious, because a nonlocal spatial operator sees beyond the local boundary stencil.

For this project:

- avoid spatial fractional Laplacians until we choose an operator definition;
- never compare spectral and integral results as if they were the same model;
- treat absorbing boundaries for nonlocal spatial operators as a research task, not a parameter tweak.

### 5.6 Stability Is Not Inherited from the Baseline

The convergence discussion in §5.3 is about accuracy; stability is a separate question and is currently unaddressed. The classical CFL condition governs the explicit leapfrog Laplacian update, but bolting a fractional-in-time memory term onto that update changes the stability picture — the memory term contributes its own constraint, and in the `1 < α < 2` wave regime the combined condition is not the bare CFL bound. This deserves its own treatment.

For this project:

- do not assume the baseline CFL number remains safe once the memory term is active;
- add an explicit stability sweep (vary `dt` at fixed `dx`) for each fractional scheme, separate from the convergence sweep;
- record the empirically stable region alongside the formal one.

## 6. The Approach I Would Take

The safest path is incremental.

**Stage 1 — Keep the fractional operator temporal.** Start with fractional time memory, not a spatial fractional Laplacian, which directly targets attenuation/dispersion while avoiding nonlocal boundary ambiguity. Use a fractional Zener-inspired law where each material owns `c_m, ρ_m, α_m, β_m, τ_σ,m, τ_ε,m` plus approximation metadata; update each cell with its own material memory law; let the ordinary Laplacian stencil do the spatial coupling.

```text
d_t^2 p = c(x)^2 Δp + memory_loss[p; material(x)] + s.
```

**Stage 2 — Implement full-history first on tiny grids.** Before optimizing, implement the honest expensive version,

```text
D_t^α u(t_n) ≈ Σ_{k=0}^n a_k u(t_{n-k}),
```

as a reference. The README already prototypes Grünwald-Letnikov weights; GL, the L1 scheme, and the L2-1σ scheme are the three standard discretizations and they differ in accuracy and in nonsmooth-data behavior, so pick one deliberately and record why. This version will be slow but tells us what any fast approximation is supposed to reproduce.

**Stage 3 — Add the SOE / diffusive approximation as an explicit approximation.** Then add `kernel(t) ≈ Σ_j w_j e^(-λ_j t)`. Do not hide it behind the same type as the full-history derivative. Store the fit interval, number of exponentials, weight signs, max kernel error, and the frequency band where the attenuation slope is acceptable.

**Stage 4 — Compare against baselines.** For every geometry, run: lossless; simple damped; Zener/SLS; constant-Q; Biot/EDFM-style if saturated sediment is involved; fractional full-history; fractional fast approximation. The fractional model earns its keep only if it predicts something the baselines do not, preferably tied to frequency dependence.

## 7. What Would Be Visually Obvious?

For a blog post, visuals matter, and some failures are visible: a circular wavefront in uniform material; a delayed/distorted wavefront through a slow block; scattering/refraction around a fast target in a slow background; damped vs lossless amplitude decay; fractional vs Zener tail behavior after the main pulse passes.

But the failures that matter most are *not* visually obvious: a wrong attenuation exponent, a phase velocity slightly off by frequency, an SOE approximation matching early time but failing late time, convergence-order collapse near `t = 0`, or boundary artifacts from a spatial fractional operator. Those need probe traces and spectra, not video.

## 8. Blog Framing

The interesting story is not "fractional equations are better." A stronger framing:

> Fractional differential equations are a compact language for memory. But numerical approximations are also memory models. Approximate carelessly and you may replace a fractional material with a finite collection of ordinary relaxation modes — without noticing.

That gives the project an arc:

1. Start with the joke premise: finding bones with sound waves.
2. Strip it to the real question: waves through heterogeneous, lossy media.
3. Show why baseline acoustics is insufficient.
4. Introduce fractional derivatives as power-law memory (and complete monotonicity as why that memory is a spectrum of relaxations).
5. Introduce spatial variation by material.
6. Hit the twist: approximations can erase the fractional behavior.
7. Build tests and visuals that expose the difference.

## 9. Practical Claims We Can Make Now

Safe to say:

- We have baseline scalar acoustic and lossy proxy models.
- We have material maps with per-cell wave speed.
- We have visual sanity scenarios.
- We have a literature-backed warning list for fractional numerics.
- We have not yet implemented a real calibrated fractional Zener solver.
- We should not claim fractional superiority until it beats Zener/constant-Q/Biot baselines on frequency-dependent probe diagnostics.

Not safe to say yet:

- that variable-order fractional models are stable in our setting;
- that SOE approximations preserve the fractional advantage;
- that short-memory truncation is acceptable;
- that spatial fractional Laplacians are straightforward on our bounded grid;
- that the reduced Biot proxy is a real poroelastic solver.

## 10. Source Review

**Fractional wave and material modeling.** Nasholm and Holm derive fractional acoustic and Zener-style wave equations from conservation laws and relaxation mechanisms, which is what keeps fractional derivatives tied to material physics rather than arbitrary curve-fitting; their comparison paper surveys the competing fractional wave equations for power-law attenuation. Baker and Banjai analyze a wave equation for lossy media obeying a frequency power law, close in spirit to what we want. Mainardi is the standard reference connecting fractional viscoelasticity, the Zener model, and the special-function/complete-monotonicity machinery; Hanyga gives the complete-monotonicity test behind the `1/2 ≤ α ≤ 1` condition for the `β(p) = C p^α` power-law attenuation example; Schilling, Song, and Vondracek is the rigorous backbone for Bernstein/completely-monotone functions. On causality, Szabo's time-domain power-law wave equation and the Waters-Mobley-Miller review of the Kramers-Kronig attenuation/dispersion link are the relevant anchors.

**Saturated sediment and experimental motivation.** Argo et al. measure sound speed in water-saturated glass beads vs porosity and frequency — directly relevant to the fish tank, and a warning that controlled granular media still show frequency dependence simple models miss. Chintada, Rau, and Goksel frame measurement spectrally (reconstructing speed of sound and attenuation rather than a single "speed"). Tsiklauri is a poroelastic/Biot-adjacent reference for fluid-saturated porous solids.

**Numerical warnings.** Jin, Lazarov, Zhou and Li, Wang, Xie give the nonsmooth-data warning (formal convergence order can fail), for subdiffusion and the fractional-wave regime respectively. Lewandowska and Kosztolowicz give the short-memory warning. Quan, Wu, Yang and Chaudhary et al. give the SOE warning (fast memory approximations must be analyzed and calibrated). Zheng gives the variable-order warning (loss of coefficient monotonicity). Lischke et al. and Borthagaray et al. give the spatial-fractional warning (boundary definitions and singularities matter).

## References

- Nasholm and Holm, "A Fractional Acoustic Wave Equation from Multiple Relaxation Loss and Conservation Laws", arXiv:1202.4251, https://arxiv.org/abs/1202.4251.
- Nasholm and Holm, "On a Fractional Zener Elastic Wave Equation", arXiv:1212.4024, https://arxiv.org/abs/1212.4024.
- Holm and Nasholm, "Comparison of fractional wave equations for power law attenuation in ultrasound and elastography", arXiv:1306.6507, https://arxiv.org/abs/1306.6507.
- Baker and Banjai, "Numerical analysis of a wave equation for lossy media obeying a frequency power law", arXiv:2012.04520, https://arxiv.org/abs/2012.04520.
- Mainardi, "Fractional Calculus and Waves in Linear Viscoelasticity: An Introduction to Mathematical Models", Imperial College Press, London, 2010 (2nd ed., World Scientific, 2022). Appendix on complete monotone and Bernstein functions.
- Hanyga, "Wave propagation in linear viscoelastic media with completely monotonic relaxation moduli", Wave Motion 50 (2013) 909-928, arXiv:1302.0402, https://arxiv.org/abs/1302.0402.
- Schilling, Song, and Vondracek, "Bernstein Functions: Theory and Applications", De Gruyter, 2nd ed., 2012.
- Szabo, "Time domain wave equations for lossy media obeying a frequency power law", J. Acoust. Soc. Am. 96 (1994) 491-500.
- Waters, Mobley, and Miller, "Causality-imposed (Kramers-Kronig) relationships between attenuation and dispersion", IEEE Trans. Ultrason. Ferroelectr. Freq. Control 52 (2005) 822-833.
- Argo IV et al., "Sound speed in water-saturated glass beads as a function of frequency and porosity", arXiv:0906.4798, https://arxiv.org/abs/0906.4798.
- Chintada, Rau, and Goksel, "Spectral Ultrasound Imaging of Speed-of-Sound and Attenuation Using an Acoustic Mirror", arXiv:2201.01435, https://arxiv.org/abs/2201.01435.
- Tsiklauri, "Phenomenological model of propagation of the elastic waves in a fluid-saturated porous solid with non-zero boundary slip velocity", arXiv:physics/0201045, https://arxiv.org/abs/physics/0201045.
- Jin, Lazarov, and Zhou, "An analysis of the L1 Scheme for the subdiffusion equation with nonsmooth data", arXiv:1501.00253, https://arxiv.org/abs/1501.00253.
- Jin, Lazarov, and Zhou, "Numerical methods for time-fractional evolution equations with nonsmooth data: a concise overview", arXiv:1805.11309, https://arxiv.org/abs/1805.11309.
- Li, Wang, and Xie, "Analysis of the L1 scheme for fractional wave equations with nonsmooth data", arXiv:1908.09145, https://arxiv.org/abs/1908.09145.
- Lewandowska and Kosztolowicz, "Numerical study of subdiffusion equation", arXiv:cond-mat/0611014, https://arxiv.org/abs/cond-mat/0611014.
- Quan, Wu, and Yang, "Long time H1-stability of fast L2-1_sigma method on general nonuniform meshes for subdiffusion equations", arXiv:2212.00453, https://arxiv.org/abs/2212.00453.
- Chaudhary, Diethelm, Farhadi, and Fuchs, "An Efficient Exponential Sum Approximation of Power-Law Kernels for Solving Fractional Differential Equation", arXiv:2508.20311, https://arxiv.org/abs/2508.20311.
- Zheng, "Numerical approximation for a nonlinear variable-order fractional differential equation via an integral equation method", arXiv:2110.04707, https://arxiv.org/abs/2110.04707.
- Lischke et al., "What Is the Fractional Laplacian?", arXiv:1801.09767, https://arxiv.org/abs/1801.09767.
- Borthagaray, Leykekhman, and Nochetto, "Local energy estimates for the fractional Laplacian", arXiv:2005.03786, https://arxiv.org/abs/2005.03786.
