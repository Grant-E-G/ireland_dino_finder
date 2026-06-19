# Fractional Numerics Pitfalls

This note tracks theory-side traps that matter for the solver. The central warning: fractional models are attractive because they compactly encode long memory and power-law loss/dispersion, and several practical approximations erase exactly that advantage.

This file is the actionable trap list. The conceptual grounding — complete monotonicity and the sum-of-exponentials connection, the Kramers-Kronig link between attenuation and dispersion, and the subdiffusion-vs-wave regime split — lives in `docs/variable_space_fractional_wave_review.md`. Here we just record what to do and what to test.

One framing fact used throughout: most of the cited numerical-analysis results were proved for *subdiffusion*, order `0 < α < 1`. This project is a *wave* problem, but the order notation splits: fractional Zener material exponents usually live in `0 < α, β <= 1`, while direct fractional-wave equations may involve effective time orders in `1 < γ < 2` or terms such as `D_t^(β + 2)`. The qualitative lessons carry over; the regularity and stability constants do not transfer automatically, and the relevant in-regime reference is usually different. Each item below is tagged with the regime its sources come from.

## 1. Short Memory Can Destroy the Model

*(sources: subdiffusion)* Fractional time derivatives are history integrals with slow power-law decay. Truncating history after a fixed window can convert a long-memory process into an ordinary finite-memory model. Lewandowska and Kosztolowicz make this point directly for subdiffusion.

Project rule:

- No fixed-window history truncation as the default.
- If truncation is added, expose it as an approximation with a tunable memory horizon and test it against full-history results.
- Any power-law attenuation claim must state the frequency/time range over which the memory approximation holds.

Local PDF: `source/pdf/lewandowska_kosztolowicz_2006_short_memory_subdiffusion.pdf`

## 2. Sum-of-Exponentials Can Become Just a Many-Relaxation Model

*(sources: general / subdiffusion stability)* Approximating a power-law kernel by a finite sum of exponentials turns expensive history convolution into recursive internal states `z_j' = -λ_j z_j + input`. That is exactly a finite collection of relaxation modes, so with too few exponentials or a poorly chosen fit interval `[δ, T]` the model is no longer broadband fractional — it is a many-relaxation Zener model. That may be a fine engineering model, but it is a different claim.

The exact kernel is completely monotone, i.e. a positive superposition of decaying exponentials (Bernstein). A faithful SOE preserves that structure: positive weights and decays keep the approximate kernel completely monotone, which is what preserves passivity. Negative weights can fit the curve while breaking the physics. Separately, Hanyga's power-law attenuation example `beta(p) = C p^alpha` gives a sharper admissibility test: complete monotonicity of the associated relaxation modulus requires roughly `1/2 ≤ α ≤ 1`. Exponents outside the relevant passive range should be flagged per material instead of treated as automatically physical.

Quan, Wu, and Yang's fast L2-1σ stability work is the warning that SOE-scheme stability is a theorem with conditions, not an automatic consequence of swapping the kernel. Chaudhary et al. frame the SOE approximation by explicit cost/error tradeoffs.

Project rule:

- Store the fitted time band `[δ, T]`, weights, decays, and kernel error with every SOE approximation.
- Test the fitted kernel against the target power law before running a wave simulation.
- Probe-spectrum check: if the SOE model no longer produces the intended attenuation slope, it has lost the fractional advantage.
- Prefer positive weights/decays; positivity ties to complete monotonicity, dissipation, and stability.
- Check that the fitted attenuation and dispersion stay Kramers-Kronig consistent over the source band (they cannot be fit independently).

Local PDFs:

- `source/pdf/quan_wu_yang_2022_fast_l2_stability_soe.pdf`
- `source/pdf/chaudhary_diethelm_farhadi_fuchs_2025_exponential_sum_power_law.pdf`

## 3. Advertised High Order Assumes Smoothness We Do Not Have

*(sources: subdiffusion — Jin/Lazarov/Zhou; wave — Li/Wang/Xie)* Fractional evolution problems have weak singular behavior near `t = 0`, so high-order schemes degrade when the solution lacks the assumed smoothness. Jin, Lazarov, and Zhou show the L1 scheme's `O(τ^(2-α))` claim depends on restrictive smoothness, with `O(τ)` for nonsmooth data; their overview collects related convolution-quadrature and L1 results. The same applies to the wave regime: Li, Wang, and Xie analyze L1 for fractional *wave* equations with nonsmooth data and propose a modified L1 scheme. Our Ricker pulses, sharp material interfaces, and abrupt source turn-on are exactly the nonsmooth cases these results target.

The discretization choice itself matters. The README prototypes Grünwald-Letnikov weights; GL, the L1 scheme, and L2-1σ differ in formal order and in nonsmooth-data robustness (GL in particular is prone to order reduction without startup correction). Pick one deliberately and record why.

Project rule:

- Do not trust nominal scheme order unless the test problem has the required regularity.
- Benchmark with both smooth manufactured signals and pulse/interface cases; estimate observed order by halving `dt` and comparing probe traces.
- Offer graded-time or startup-correction options before claiming high order.
- Keep the fast smoke tests, but add a separate convergence suite for the fractional-wave variants with discontinuous material maps and compact pulses.
- Compare full-history, SOE, and any short-memory approximations on the same probe traces.

Local PDFs:

- `source/pdf/jin_lazarov_zhou_2015_l1_nonsmooth_data.pdf`
- `source/pdf/jin_lazarov_zhou_2018_nonsmooth_data_overview.pdf`
- `source/pdf/li_wang_xie_2019_l1_fractional_wave_nonsmooth.pdf`

## 4. Variable Order Breaks Some Standard Proof Machinery

*(sources: general variable-order)* Space-varying fractional order by material is not constant-order with `α` replaced by `α(x)`. Zheng's variable-order work shows standard schemes can lose monotonicity of discretization coefficients, so existing analysis does not apply directly — a real warning for stability and convergence if `α` jumps across interfaces.

Project rule:

- Start with piecewise-constant orders per material and fixed interface tests.
- Do not interpolate fractional order smoothly without a physical reason and a test.
- Add a two-material interface benchmark before any 2D/3D production visuals.
- Treat coefficient monotonicity/positivity as a diagnostic, not an implementation detail.

Local PDF: `source/pdf/zheng_2021_variable_order_integral_equation.pdf`

## 5. Spatial Fractional Operators Have Boundary-Definition Traps

*(sources: spatial fractional Laplacian)* If we later add a spatial fractional Laplacian, boundary conditions become part of the operator definition. Lischke et al. show Riesz/integral, spectral, directional, and horizon-based fractional Laplacians are not interchangeable on bounded domains — an integral definition can need exterior values, a spectral one uses local boundary data. Borthagaray, Leykekhman, and Nochetto show boundary-singularity effects for the integral version that degrade global convergence. Two implementations can disagree near absorbing boundaries even though both are called "the fractional Laplacian."

Project rule:

- Choose and document the operator definition before implementing spatial fractional operators.
- Do not mix spectral and integral results in the same validation table.
- Treat absorbing boundaries/PMLs as an open research item for nonlocal spatial operators.
- Add interior-vs-boundary error diagnostics if spatial fractional operators enter the solver.

Local PDFs:

- `source/pdf/lischke_pang_gulian_2018_what_is_fractional_laplacian.pdf`
- `source/pdf/borthagaray_leykekhman_nochetto_2020_fractional_laplacian_boundary_singularity.pdf`

## 6. Stability Is Not Inherited from the Baseline CFL

*(sources: general)* The items above are mostly about accuracy. Stability is separate and currently unverified. The classical CFL condition governs the explicit leapfrog Laplacian update, but adding a fractional-in-time memory term changes the stability picture: the memory term contributes its own constraint, and in the `1 < α < 2` regime the combined bound is not the bare CFL number. Do not assume the baseline's safe Courant number stays safe once memory is active.

Project rule:

- Run a dedicated stability sweep (vary `dt` at fixed `dx`) for each fractional scheme, separate from the convergence sweep.
- Record the empirically stable region alongside any formal bound.
- Re-check stability whenever the SOE weights, memory horizon, or material order set changes.

## Implementation Checklist

Before merging a real fractional solver:

- [ ] Define the derivative/operator precisely: Caputo, Riemann-Liouville, fractional Zener, spatial fractional Laplacian, etc.
- [ ] State the discretization (GL, L1, L2-1σ, convolution quadrature) and why.
- [ ] State whether the model is full-history, SOE, diffusive representation, short-memory, or another approximation.
- [ ] Store approximation metadata: order, time band, memory length, SOE weights/decays, and kernel error.
- [ ] Kernel-fit test against the target power law.
- [ ] Convergence test with smooth manufactured data.
- [ ] Convergence test with a nonsmooth pulse/interface case.
- [ ] Stability sweep with the memory term active (not just the baseline CFL check).
- [ ] Probe-spectrum attenuation-slope test, with a Kramers-Kronig consistency check on attenuation vs dispersion.
- [ ] Per-material flag for exponent choices that fall outside the passive range for the selected model.
- [ ] Visual scenario where the fractional model must visibly differ from damped/Zener baselines.
- [ ] Document when the approximation degenerates into a finite relaxation model.

## References

- Lewandowska and Kosztolowicz, "Numerical study of subdiffusion equation", arXiv:cond-mat/0611014, https://arxiv.org/abs/cond-mat/0611014.
- Jin, Lazarov, and Zhou, "An analysis of the L1 Scheme for the subdiffusion equation with nonsmooth data", arXiv:1501.00253, https://arxiv.org/abs/1501.00253.
- Jin, Lazarov, and Zhou, "Numerical methods for time-fractional evolution equations with nonsmooth data: a concise overview", arXiv:1805.11309, https://arxiv.org/abs/1805.11309.
- Li, Wang, and Xie, "Analysis of the L1 scheme for fractional wave equations with nonsmooth data", arXiv:1908.09145, https://arxiv.org/abs/1908.09145.
- Zheng, "Numerical approximation for a nonlinear variable-order fractional differential equation via an integral equation method", arXiv:2110.04707, https://arxiv.org/abs/2110.04707.
- Quan, Wu, and Yang, "Long time H1-stability of fast L2-1_sigma method on general nonuniform meshes for subdiffusion equations", arXiv:2212.00453, https://arxiv.org/abs/2212.00453.
- Chaudhary, Diethelm, Farhadi, and Fuchs, "An Efficient Exponential Sum Approximation of Power-Law Kernels for Solving Fractional Differential Equation", arXiv:2508.20311, https://arxiv.org/abs/2508.20311.
- Lischke et al., "What Is the Fractional Laplacian?", arXiv:1801.09767, https://arxiv.org/abs/1801.09767.
- Borthagaray, Leykekhman, and Nochetto, "Local energy estimates for the fractional Laplacian", arXiv:2005.03786, https://arxiv.org/abs/2005.03786.
- Hanyga, "Wave propagation in linear viscoelastic media with completely monotonic relaxation moduli", Wave Motion 50 (2013) 909-928, arXiv:1302.0402, https://arxiv.org/abs/1302.0402.
- Mainardi, "Fractional Calculus and Waves in Linear Viscoelasticity: An Introduction to Mathematical Models", Imperial College Press, London, 2010 (2nd ed., World Scientific, 2022).
- Waters, Mobley, and Miller, "Causality-imposed (Kramers-Kronig) relationships between attenuation and dispersion", IEEE Trans. Ultrason. Ferroelectr. Freq. Control 52 (2005) 822-833.
