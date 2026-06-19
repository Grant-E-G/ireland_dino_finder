# Variable-Space Fractional Wave Models: Review Notes

These are working notes for a longer blog post. They assume comfort with analysis, dynamical systems, operators, and asymptotics, but not much recent acoustics or applied numerical PDE background.

The motivating question is:

> Can a wave simulation use spatially varying fractional differential equations to model how sound behaves in different materials, without throwing away the very memory/dispersion effects that made fractional equations interesting?

The short answer from the sources is: maybe, but only if we are careful. Fractional models are not magic damping terms. They are nonlocal-in-time or nonlocal-in-space operators with real numerical and modeling traps. Some common approximations can collapse the model back into a finite-memory relaxation system, which may be useful, but is no longer the thing we thought we were testing.

## 1. Physics Setup: What Is Being Modeled?

For a scalar acoustic approximation, the clean baseline is pressure `p(x, t)` satisfying something like

```text
d_t^2 p = c(x)^2 Δp + s.
```

Here:

- `p` is pressure or pressure-like wave amplitude;
- `c(x)` is local wave speed;
- `s` is a source term;
- `Δ` is the spatial Laplacian.

This equation is conservative in the ideal homogeneous case. A pulse travels, reflects, refracts, and interferes, but it does not intrinsically lose energy to heat, pore-fluid motion, microstructure, or scattering below the grid scale.

Real materials are not that clean. Water, sand, saturated sediment, rock, and bone-like materials differ not only by wave speed, but by attenuation and dispersion:

- attenuation: amplitude decreases as the wave propagates;
- dispersion: phase velocity depends on frequency;
- scattering: energy moves into unresolved directions/modes;
- poroelastic loss: pore fluid and solid frame move relative to one another;
- relaxation: internal degrees of freedom lag behind stress/strain.

Experimentally, attenuation is often fit by a power law over a frequency band:

```text
attenuation(f) = a0 f^y,
```

where `y` is not generally an integer. This is the door through which fractional equations enter.

## 2. Why Fractional Derivatives Appear

The usual integer-order models give a small menu of frequency behaviors. A single viscous damping term such as

```text
d_t^2 p + γ d_t p = c^2 Δp
```

does not give arbitrary broadband power-law attenuation and dispersion. A single Zener or standard-linear-solid relaxation gives one characteristic time scale. Several relaxations give several time scales.

A fractional derivative is a way of encoding a continuum-like distribution of relaxation times. For a Caputo derivative of order `0 < α < 1`,

```text
D_t^α u(t) = 1 / Γ(1 - α) ∫_0^t (t - τ)^(-α) u'(τ) dτ.
```

The important feature is the kernel:

```text
(t - τ)^(-α).
```

That kernel decays slowly. The present state depends on the whole past, with a power-law memory. In frequency space, fractional derivatives behave like powers of `iω`:

```text
D_t^α  ↔  (iω)^α.
```

That is the compact connection to non-integer power-law attenuation and dispersion.

For materials, the more physical story is often phrased through a fractional Zener constitutive relation:

```text
σ + τ_σ^α D_t^α σ
  = E0 (ε + τ_ε^β D_t^β ε).
```

This says stress and strain are related through a memory law. The material is neither purely elastic nor simply viscous. The fractional exponents tune the way the material remembers deformation.

## 3. The Space-Varying Idea

The project idea is not only "use a fractional derivative." It is:

```text
α = α(x)
β = β(x)
τ_σ = τ_σ(x)
τ_ε = τ_ε(x)
c = c(x)
ρ = ρ(x)
```

Each material gets its own fractional behavior. Water might be close to ordinary acoustics. Saturated sand might have strong poroelastic loss. Bone-like material might have a different attenuation exponent and velocity. A grid cell's material id determines the local parameters.

The conservative first implementation choice is piecewise constant material parameters:

```text
α(x) = α_j  when x is in material region j.
```

This matches how we already represent material maps. It also avoids pretending we know how to physically interpolate fractional order across a sharp boundary.

The first variable-space fractional wave model should therefore look more like a material-indexed operator family than a globally smooth coefficient field:

```text
cell i has material m(i)
cell i uses α_m, τ_m, c_m, attenuation parameters, history approximation m.
```

That has a dynamical-systems flavor: each material region carries local internal memory state, and the pressure field couples those local memory systems through the spatial stencil.

## 4. The Baselines We Need Before Fractional Models

The baselines are not busywork. They are how we keep ourselves honest.

### 4.1 Lossless Acoustic Wave

```text
d_t^2 p = c(x)^2 Δp + s.
```

This tests geometry, wave speed, finite-difference stability, and visual sanity. If this fails, nothing fractional matters.

What it can do:

- arrival-time estimates;
- reflection/refraction from wave-speed contrasts;
- basic target visibility;
- code-level debugging.

What it cannot do:

- intrinsic attenuation;
- frequency-dependent loss;
- frequency-dependent phase velocity;
- unresolved scattering or poroelastic effects.

### 4.2 Simple Damped Acoustic Wave

```text
d_t^2 p + γ(x) d_t p = c(x)^2 Δp + s.
```

This is the minimal lossy model. It is useful for checking that energy decays, but it is not a serious broadband attenuation law.

The important comparison is: if the fractional model does not visibly or spectrally differ from this, the fractional machinery may be dead weight.

### 4.3 Standard Linear Solid / Zener

A single relaxation mechanism can be expressed through a frequency-dependent modulus:

```text
M(ω) = M_∞ (1 + iωτ_ε) / (1 + iωτ_σ).
```

This has physical meaning: there is a characteristic relaxation time. It gives causal attenuation and dispersion, but only over a limited shape. Multiple Zener mechanisms can approximate broader behavior.

This baseline matters because many numerical fractional approximations effectively turn into "some number of Zener mechanisms." That is not automatically bad, but it changes the interpretation.

### 4.4 Constant-Q / Kjartansson-Style

Seismology often uses quality factor `Q`, roughly measuring energy lost per cycle:

```text
Q^-1 ≈ energy lost per cycle / (2π stored energy).
```

Constant-Q models are appealing because they encode broadband attenuation in a compact way. But causal constant-Q behavior is subtle. Low- and high-frequency cutoffs matter.

This is a useful baseline because it sits between "simple damping" and "fractional constitutive law."

### 4.5 Biot / Poroelastic Models

For saturated sand or sediment, scalar acoustics is physically too thin. Biot theory couples solid-frame displacement and pore-fluid displacement:

```text
solid momentum: ρ_11 d_t^2 u + ρ_12 d_t^2 U = div(σ)
fluid momentum: ρ_12 d_t^2 u + ρ_22 d_t^2 U + b d_t(U - u) = -∇p_f.
```

This is directly relevant to the fish-tank idea: sand plus water plus a bone-like inclusion.

The cost is parameter burden. Permeability, tortuosity, porosity, frame moduli, and pore geometry are not decoration; they drive the model. A bad Biot parameter set can be less honest than a humble scalar baseline.

Argo et al.'s water-saturated glass-bead measurements are useful here. They report porosity-dependent sound speed and negative dispersion above about 550 kHz. Their result is a warning: even Biot-derived effective-density-fluid models may match porosity trends while missing frequency trends.

## 5. What the Fractional Literature Warns Us About

The strongest theme from the numerical-analysis sources is that fractional equations are fragile under approximation. The model's advantage comes from long memory and nonlocality. Practical approximations often attack exactly those properties.

### 5.1 Short Memory Can Kill the Point

A fractional derivative remembers the whole past with a power-law tail. If we truncate history to a fixed recent window, we may no longer have a fractional model in the intended sense.

Lewandowska and Kosztolowicz explicitly warn against using the short-memory principle for subdiffusion. The details are diffusion rather than wave propagation, but the modeling lesson transfers: a long-memory process is not faithfully represented by an arbitrary finite memory window.

For this project:

- no fixed short-memory truncation as the default;
- if used, it must be labeled and tested against full-history behavior;
- visual similarity is not enough; compare probe spectra and attenuation slopes.

### 5.2 Sum-of-Exponentials Is Useful but Dangerous

A common computational trick is to approximate the power-law kernel by

```text
t^(-α) ≈ Σ_j w_j exp(-λ_j t).
```

This turns the history integral into recursive internal states:

```text
z_j'(t) = -λ_j z_j(t) + input(t).
```

That is computationally attractive. It is also exactly a finite collection of relaxation modes.

So the danger is conceptual: with too few exponentials, or with a poorly chosen fit interval, we no longer have a broadband fractional memory law. We have a many-relaxation Zener-like model.

That might be a perfectly good engineering model. It is just not the same claim.

For this project:

- store the approximation interval `[δ, T]`;
- store weights, decays, and kernel-fit error;
- require positivity/stability checks when the theory needs them;
- test whether the resulting attenuation still follows the intended power law over the source band.

### 5.3 High-Order Claims Often Require Smoothness We Do Not Have

Fractional evolution equations often have weak singularities near `t = 0`. Even with smooth forcing, the solution may not be as regular as an integer-order PDE solution. Jin, Lazarov, and Zhou show that classical convergence claims for the L1 scheme depend on restrictive smoothness; nonsmooth data can degrade the observed rate.

This matters because our simulations intentionally include:

- Ricker pulses;
- sharp material interfaces;
- abrupt material parameter jumps;
- potentially non-smooth source turn-on.

So a scheme advertised as `O(τ^(2-α))` may behave closer to first order in the cases we actually care about.

For this project:

- convergence tests must include nonsmooth pulse/interface cases;
- smooth manufactured tests are necessary but insufficient;
- startup corrections or graded time meshes may be needed;
- any blog claim about accuracy should distinguish formal order from observed order.

### 5.4 Variable Order Is Not a Text Substitution

Replacing `α` by `α(x)` looks harmless:

```text
D_t^α  →  D_t^{α(x)}.
```

But variable-order fractional equations can break proof structures that work for constant order. Zheng's variable-order work notes that discretization coefficients can lose monotonicity, which blocks standard numerical analysis techniques.

This is directly relevant to material maps. If `α` jumps at a boundary, the memory law changes discontinuously across neighboring cells.

For this project:

- begin with piecewise-constant material orders;
- treat monotonicity/positivity of discrete weights as a diagnostic;
- add a two-material interface benchmark before doing elaborate visuals;
- avoid arbitrary smoothing of `α(x)` unless it has physical or numerical justification.

### 5.5 Spatial Fractional Operators Are Boundary-Condition Minefields

So far the project mainly points toward fractional time derivatives. But "space fractional" can also mean fractional Laplacians:

```text
(-Δ)^s u.
```

On all of `R^n`, this is already nonlocal. On a bounded simulation box, the phrase "the fractional Laplacian" is not enough. Lischke et al. explain that integral, spectral, directional, and other definitions differ. Boundary behavior is part of the operator. Borthagaray, Leykekhman, and Nochetto show boundary singularity effects for the integral fractional Laplacian.

For wave simulation with absorbing boundaries, this is a serious issue. A nonlocal spatial operator sees beyond the local boundary stencil.

For this project:

- avoid spatial fractional Laplacians until we choose an operator definition;
- do not compare spectral and integral fractional Laplacian results as if they were the same model;
- treat absorbing boundaries for nonlocal spatial operators as a research task, not a parameter tweak.

## 6. The Approach I Would Take

The safest path is incremental.

### Stage 1: Keep the Fractional Operator Temporal

Start with fractional time memory, not a spatial fractional Laplacian. This directly targets attenuation and dispersion while avoiding nonlocal boundary-condition ambiguity.

The model family is:

```text
d_t^2 p = c(x)^2 Δp + memory_loss[p; material(x)] + s.
```

More concretely, use a fractional Zener-inspired law where each material owns:

```text
material m:
  c_m
  ρ_m
  α_m
  β_m
  τ_σ,m
  τ_ε,m
  approximation metadata
```

At each cell:

```text
update pressure using that cell's material memory law.
```

Spatial coupling still comes from the ordinary Laplacian stencil.

### Stage 2: Implement Full-History First on Tiny Grids

Before optimizing, implement the honest expensive version on very small grids:

```text
D_t^α u(t_n) ≈ Σ_{k=0}^n a_k u(t_{n-k})
```

This gives a reference. It will be slow, but it will tell us what the approximation is supposed to reproduce.

### Stage 3: Add SOE / Diffusive Approximation as an Explicit Approximation

Then add the fast representation:

```text
kernel(t) ≈ Σ_j w_j e^{-λ_j t}.
```

Do not hide this behind the same type as the full-history derivative. Store:

- fit interval;
- number of exponentials;
- weight signs;
- max kernel error;
- frequency band where attenuation slope is acceptable.

### Stage 4: Compare Against Baselines

For every geometry, run:

1. lossless acoustic;
2. simple damped acoustic;
3. Zener/SLS;
4. constant-Q;
5. Biot/EDFM-style if saturated sediment is involved;
6. fractional full-history;
7. fractional fast approximation.

The fractional model earns its keep only if it predicts something the baselines do not, preferably something tied to frequency dependence.

## 7. What Would Be Visually Obvious?

For a blog post, visuals matter. Some failures should be visible.

Good visual sanity checks:

- circular wavefront in uniform material;
- delayed/distorted wavefront through a slow block;
- scattering/refraction around a fast circular target in a slow background;
- damped vs lossless amplitude decay;
- fractional vs Zener tail behavior after the main pulse passes.

But some important failures will not be visually obvious:

- wrong attenuation exponent;
- phase velocity slightly wrong by frequency;
- SOE approximation matching early time but failing late time;
- convergence order collapse near `t = 0`;
- boundary artifacts from a spatial fractional operator.

Those need probe traces and spectra.

## 8. Blog Framing

The interesting story is not "fractional equations are better." A stronger framing is:

> Fractional differential equations are a compact language for memory. But numerical approximations are also memory models. If we approximate carelessly, we may accidentally replace a fractional material with a finite collection of ordinary relaxation modes.

That gives the project a useful arc:

1. Start with a joke premise: finding bones with sound waves.
2. Strip it down to the real question: waves through heterogeneous, lossy media.
3. Show why baseline acoustics is insufficient.
4. Introduce fractional derivatives as power-law memory.
5. Introduce spatial variation by material.
6. Hit the twist: approximations can erase the fractional behavior.
7. Build tests and visuals that expose the difference.

## 9. Practical Claims We Can Make Now

It is safe to say:

- We have baseline scalar acoustic and lossy proxy models.
- We have material maps with per-cell wave speed.
- We have visual sanity scenarios.
- We have a literature-backed warning list for fractional numerics.
- We have not yet implemented a real calibrated fractional Zener solver.
- We should not claim fractional superiority until it beats Zener/constant-Q/Biot baselines on frequency-dependent probe diagnostics.

It is not safe to say yet:

- that variable-order fractional models are stable in our setting;
- that SOE approximations preserve the fractional advantage;
- that short-memory truncation is acceptable;
- that spatial fractional Laplacians are straightforward on our bounded grid;
- that the reduced Biot proxy is a real poroelastic solver.

## 10. Source Review

The sources split into three groups.

### Fractional Wave and Material Modeling

Nasholm and Holm derive fractional acoustic and Zener-style wave equations from conservation laws and relaxation mechanisms. Their work is useful because it connects fractional derivatives to material relaxation rather than treating them as arbitrary curve-fitting. Holm and Nasholm compare fractional wave equations for power-law attenuation in ultrasound and elastography.

Baker and Banjai analyze a wave equation for lossy media obeying a frequency power law. This is close in spirit to what we want: a wave model where the loss law is tied to measured frequency behavior.

### Saturated Sediment and Experimental Motivation

Argo et al. measure sound speed in water-saturated glass beads as a function of porosity and frequency. This is highly relevant to the fish-tank experiment. Their data show that even controlled granular materials can have frequency-dependent behavior that simple models do not fully capture.

Chintada, Rau, and Goksel are useful for measurement framing: reconstructing speed of sound and attenuation spectrally, rather than treating "the speed" as a single number.

Tsiklauri provides a poroelastic/Biot-adjacent reference point for fluid-saturated porous solids.

### Numerical Warnings

Jin, Lazarov, Zhou, Li, Wang, and Xie provide the nonsmooth-data warning: formal convergence order can fail for realistic fractional evolution problems.

Lewandowska and Kosztolowicz provide the short-memory warning: truncating memory can destroy the model class.

Quan, Wu, Yang, and Chaudhary et al. provide the SOE warning: fast memory approximations are powerful but must be analyzed and calibrated.

Zheng provides the variable-order warning: changing `α` across space can break coefficient monotonicity and standard analysis.

Lischke et al. and Borthagaray et al. provide the spatial fractional warning: boundary definitions and singularities matter.

## References

- Nasholm and Holm, "A Fractional Acoustic Wave Equation from Multiple Relaxation Loss and Conservation Laws", arXiv:1202.4251, https://arxiv.org/abs/1202.4251.
- Nasholm and Holm, "On a Fractional Zener Elastic Wave Equation", arXiv:1212.4024, https://arxiv.org/abs/1212.4024.
- Holm and Nasholm, "Comparison of fractional wave equations for power law attenuation in ultrasound and elastography", arXiv:1306.6507, https://arxiv.org/abs/1306.6507.
- Baker and Banjai, "Numerical analysis of a wave equation for lossy media obeying a frequency power law", arXiv:2012.04520, https://arxiv.org/abs/2012.04520.
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
