# Fractional Numerics Pitfalls

This note tracks theory-side traps that matter for this project. The central warning is simple: fractional models are attractive because they compactly encode long memory and power-law loss/dispersion. Several practical approximations can erase exactly that advantage.

## 1. Short Memory Can Destroy the Model

Fractional time derivatives are history integrals with slow power-law decay. If we simply truncate history after a fixed window, we may convert a long-memory process into an ordinary finite-memory model. Lewandowska and Kosztolowicz make this point directly for subdiffusion: the process is long-memory, and the short-memory principle should not be used for that case.

Project rule:

- Do not use fixed-window history truncation as the default fractional derivative.
- If we add truncation, expose it as an approximation with a tunable memory horizon and test it against full-history results.
- Any claim about power-law attenuation must include a frequency/time range over which the memory approximation is valid.

Local PDF:

- `source/pdf/lewandowska_kosztolowicz_2006_short_memory_subdiffusion.pdf`

## 2. Sum-of-Exponentials Can Become Just a Many-Relaxation Model

Approximating a power-law kernel by a finite sum of exponentials is useful because it turns expensive history convolution into a set of recursive internal states. The catch is conceptual and numerical:

- too few exponentials give a kernel that behaves like a small number of Debye/Zener relaxations, not a broadband fractional power law;
- the approximation is valid only on a chosen time interval `[delta, T]`;
- kernel compression can change the effective memory law;
- stability proofs may depend on positivity/semidefinite properties of the approximated convolution weights.

This is not a reason to avoid sum-of-exponentials. It is a reason to treat it as a calibrated approximation, not as "the" fractional derivative. Quan, Wu, and Yang's fast L2-1_sigma stability work is a useful warning: stability of fast SOE schemes is a theorem with conditions, not an automatic consequence of replacing the kernel. Chaudhary, Diethelm, Farhadi, and Fuchs are useful for implementation because they frame SOE approximation by explicit cost/error tradeoffs.

Project rule:

- Store the fitted time band and kernel error with every SOE approximation.
- Test the fitted kernel against the target power law before running a wave simulation.
- Include a probe-spectrum check: if the SOE model no longer produces the intended attenuation slope, it has lost the fractional advantage.
- Prefer positive weights/decays where possible, because positivity is tied to dissipation and stability arguments.

Local PDFs:

- `source/pdf/quan_wu_yang_2022_fast_l2_stability_soe.pdf`
- `source/pdf/chaudhary_diethelm_farhadi_fuchs_2025_exponential_sum_power_law.pdf`

## 3. Advertised High Order Often Assumes Smooth Solutions We Do Not Have

Fractional evolution problems commonly have weak singular behavior near `t = 0`. High-order schemes can degrade sharply when the solution lacks the smoothness assumed by the analysis. Jin, Lazarov, and Zhou show that the popular L1 scheme's classical `O(tau^(2-alpha))` claim depends on restrictive smoothness; for nonsmooth data they establish `O(tau)` behavior. Their overview paper collects related results for convolution quadrature, L1-type schemes, and nonsmooth data.

This matters for wave simulations because Ricker pulses, material interfaces, and abrupt initial/source conditions are not globally smooth in the sense required by many clean convergence claims.

Project rule:

- Do not trust nominal scheme order unless the test problem has the required regularity.
- Include graded-time or startup-correction options before claiming high order.
- Benchmark with both smooth manufactured signals and pulse/interface cases.
- Track convergence empirically: halve `dt`, compare probe traces, and estimate observed order.

Local PDFs:

- `source/pdf/jin_lazarov_zhou_2015_l1_nonsmooth_data.pdf`
- `source/pdf/jin_lazarov_zhou_2018_nonsmooth_data_overview.pdf`
- `source/pdf/li_wang_xie_2019_l1_fractional_wave_nonsmooth.pdf`

## 4. Variable Order Breaks Some Standard Proof Machinery

This project is specifically interested in space-varying fractional order by material. Variable-order fractional derivatives are not just constant-order derivatives with `alpha` replaced by `alpha(x)`.

Zheng's variable-order work points out that standard approximation schemes can lose monotonicity of discretization coefficients, so existing numerical analysis techniques do not apply directly. That is a major warning for stability and convergence if we let `alpha` jump across material interfaces.

Project rule:

- Start with piecewise-constant orders per material and fixed interface tests.
- Do not interpolate fractional order smoothly unless there is a physical reason and a test.
- Add a two-material interface benchmark before any 2D/3D production visuals.
- Treat coefficient monotonicity/positivity as a diagnostic, not an implementation detail.

Local PDF:

- `source/pdf/zheng_2021_variable_order_integral_equation.pdf`

## 5. Spatial Fractional Operators Have Boundary-Definition Traps

If we later add a spatial fractional Laplacian, boundary conditions become part of the operator definition. Lischke et al. show that Riesz/integral, spectral, directional, and horizon-based fractional Laplacians are not interchangeable on bounded domains. For example, an integral definition can require values on the exterior of the domain, while a spectral definition uses standard local boundary data. Borthagaray, Leykekhman, and Nochetto also show boundary singularity effects for the integral fractional Laplacian that degrade global convergence.

This matters for wave visuals: two implementations can look different near absorbing boundaries even if both are called "fractional Laplacian."

Project rule:

- Before implementing spatial fractional operators, choose and document the operator definition.
- Do not mix spectral and integral fractional Laplacian results in the same validation table.
- Treat absorbing boundaries/PMLs as an open research item for nonlocal spatial operators.
- Add interior-vs-boundary error diagnostics if spatial fractional operators enter the solver.

Local PDFs:

- `source/pdf/lischke_pang_gulian_2018_what_is_fractional_laplacian.pdf`
- `source/pdf/borthagaray_leykekhman_nochetto_2020_fractional_laplacian_boundary_singularity.pdf`

## 6. Fractional Wave Equations Need Nonsmooth-Data Tests Too

Fractional diffusion results are not automatically wave results, but the same warning applies: source pulses and material jumps can reduce observed convergence. Li, Wang, and Xie analyze L1 schemes for fractional wave equations with nonsmooth data and propose a modified L1 scheme with better behavior.

Project rule:

- Keep the current fast smoke tests, but add a separate convergence suite for fractional wave variants.
- Include discontinuous material maps and compact pulses, not only smooth manufactured solutions.
- Compare full-history, SOE, and any short-memory approximations on the same probe traces.

Local PDF:

- `source/pdf/li_wang_xie_2019_l1_fractional_wave_nonsmooth.pdf`

## Implementation Checklist

Before merging a real fractional solver:

- [ ] Define the derivative/operator precisely: Caputo, Riemann-Liouville, fractional Zener, spatial fractional Laplacian, etc.
- [ ] State whether the model is full-history, SOE, diffusive representation, short-memory, or another approximation.
- [ ] Store approximation metadata: order, time band, memory length, SOE weights/decays, and kernel error.
- [ ] Add a kernel-fit test against the target power law.
- [ ] Add a convergence test with smooth manufactured data.
- [ ] Add a nonsmooth pulse/interface test.
- [ ] Add a probe-spectrum attenuation-slope test.
- [ ] Add a visual scenario where the fractional model must visibly differ from damped/Zener baselines.
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
