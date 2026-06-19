#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
pub enum WaveModel {
    /// README model 1: d_t^2 p = c^2 laplacian(p) + s.
    LosslessAcoustic,
    /// README model 2: d_t^2 p + gamma d_t p = c^2 laplacian(p) + s.
    LinearDampedAcoustic { gamma: f32 },
    /// README model 3, reduced SLS/Zener proxy with one relaxed strain/Laplacian memory.
    StandardLinearSolid {
        damping_gamma: f32,
        relaxation_time_s: f32,
        relaxation_strength: f32,
    },
    /// README model 4, band-limited constant-Q proxy around a reference frequency.
    ConstantQ {
        q: f32,
        reference_freq_hz: f32,
        dispersion_strength: f32,
    },
    /// README model 5, reduced poroelastic proxy with pore-drag memory.
    ReducedBiotPoroelastic {
        drag_gamma: f32,
        relaxation_time_s: f32,
        pore_coupling: f32,
    },
}

pub fn update_linear_damped(
    p: f32,
    p_prev: f32,
    lap: f32,
    coef_x: f32,
    dt: f32,
    gamma: f32,
) -> f32 {
    let half_step_damping = 0.5 * gamma * dt;
    (2.0 * p - (1.0 - half_step_damping) * p_prev + coef_x * lap) / (1.0 + half_step_damping)
}

pub fn update_pressure(
    model: WaveModel,
    p: f32,
    p_prev: f32,
    lap: f32,
    aux: f32,
    coef_x: f32,
    dt: f32,
) -> (f32, f32) {
    let velocity = (p - p_prev) / dt;

    match model {
        WaveModel::LosslessAcoustic => (2.0 * p - p_prev + coef_x * lap, aux),
        WaveModel::LinearDampedAcoustic { gamma } => {
            (update_linear_damped(p, p_prev, lap, coef_x, dt, gamma), aux)
        }
        WaveModel::StandardLinearSolid {
            damping_gamma,
            relaxation_time_s,
            relaxation_strength,
        } => {
            assert!(
                (0.0..=1.0).contains(&relaxation_strength),
                "relaxation_strength must be in [0, 1]"
            );
            let relaxed_lap = relax_toward(aux, lap, dt, relaxation_time_s);
            let effective_lap =
                (1.0 - relaxation_strength) * lap + relaxation_strength * relaxed_lap;
            let p_next =
                2.0 * p - p_prev + coef_x * effective_lap - dt * dt * damping_gamma * velocity;
            (p_next, relaxed_lap)
        }
        WaveModel::ConstantQ {
            q,
            reference_freq_hz,
            dispersion_strength,
        } => {
            assert!(q > 0.0, "constant-Q model requires q > 0");
            assert!(
                reference_freq_hz > 0.0,
                "constant-Q model requires reference_freq_hz > 0"
            );
            assert!(
                (0.0..=1.0).contains(&dispersion_strength),
                "dispersion_strength must be in [0, 1]"
            );
            let gamma = 2.0 * std::f32::consts::PI * reference_freq_hz / q;
            let tau = 1.0 / (2.0 * std::f32::consts::PI * reference_freq_hz);
            let relaxed_lap = relax_toward(aux, lap, dt, tau);
            let effective_lap = lap + dispersion_strength * (relaxed_lap - lap);
            (
                update_linear_damped(p, p_prev, effective_lap, coef_x, dt, gamma),
                relaxed_lap,
            )
        }
        WaveModel::ReducedBiotPoroelastic {
            drag_gamma,
            relaxation_time_s,
            pore_coupling,
        } => {
            let pore_velocity = relax_toward(aux, pore_coupling * velocity, dt, relaxation_time_s);
            let relative_velocity = velocity - pore_velocity;
            let p_next = 2.0 * p - p_prev + coef_x * lap - dt * dt * drag_gamma * relative_velocity;
            (p_next, pore_velocity)
        }
    }
}

pub fn relax_toward(current: f32, target: f32, dt: f32, relaxation_time_s: f32) -> f32 {
    assert!(
        relaxation_time_s > 0.0,
        "relaxation_time_s must be positive"
    );
    let blend = (dt / relaxation_time_s).clamp(0.0, 1.0);
    current + blend * (target - current)
}

/// Prototype Grünwald-Letnikov coefficients for a fractional derivative.
///
/// These are not yet wired into the wave solver. They give us a tested place
/// to start before deciding how much history to keep per material region.
#[allow(dead_code)]
pub fn grunwald_letnikov_weights(order: f32, n_terms: usize) -> Vec<f32> {
    assert!(order > 0.0, "fractional order must be positive");
    let mut weights = Vec::with_capacity(n_terms);
    if n_terms == 0 {
        return weights;
    }

    weights.push(1.0);
    for k in 1..n_terms {
        let previous = weights[k - 1];
        weights.push(-previous * (order - (k as f32 - 1.0)) / k as f32);
    }
    weights
}
