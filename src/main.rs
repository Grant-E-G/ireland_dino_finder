#[allow(dead_code)]
mod hdf5_helper;

use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};

fn idx(x: usize, z: usize, nx: usize) -> usize {
    z * nx + x
}

#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
enum WaveModel {
    /// README model 1: d_t^2 p = c^2 laplacian(p) + s.
    LosslessAcoustic,
    /// README model 2: d_t^2 p + gamma d_t p = c^2 laplacian(p) + s.
    LinearDampedAcoustic { gamma: f32 },
}

#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
enum Source {
    Ricker { freq_hz: f32 },
    Impulse { step: usize, amplitude: f32 },
}

struct SimParams {
    nx: usize,
    nz: usize,
    dx: f32,
    dz: f32,
    dt: f32,
    n_steps: usize,
    c0: f32,
    source_x: usize,
    source_z: usize,
    source: Source,
    absorb_width: usize,
    absorb_max: f32,
    model: WaveModel,
}

impl SimParams {
    fn default() -> Self {
        // Units are arbitrary here; think "meters" and "seconds".
        let nx = 200;
        let nz = 200;
        let dx = 1.0;
        let dz = 1.0;
        let c0 = 2000.0; // m/s
        let dt = 0.00025; // s

        // Courant number r = c dt / dx; for 2D with 5-point Laplacian, r < ~1 / sqrt(2)
        let r = c0 * dt / dx;
        println!("Courant number r = {}", r);

        Self {
            nx,
            nz,
            dx,
            dz,
            dt,
            n_steps: 2000,
            c0,
            source_x: nx / 2,
            source_z: nz / 2,
            source: Source::Ricker { freq_hz: 100.0 }, // Hz (very low; just for demo)
            absorb_width: 20,
            absorb_max: 0.015,
            model: WaveModel::LosslessAcoustic,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Probe {
    x: usize,
    z: usize,
}

struct SimOutput {
    final_pressure: Vec<f32>,
    #[allow(dead_code)]
    probe_traces: Vec<Vec<f32>>,
}

struct Fields {
    p: Vec<f32>,
    p_prev: Vec<f32>,
    p_next: Vec<f32>,
    damping: Vec<f32>,
}

impl Fields {
    fn new(nx: usize, nz: usize) -> Self {
        let n = nx * nz;
        Self {
            p: vec![0.0; n],
            p_prev: vec![0.0; n],
            p_next: vec![0.0; n],
            damping: vec![1.0; n],
        }
    }
}

/// Precompute a simple sponge-like damping coefficient near boundaries.
///
/// damping(x,z) in (0,1]; 1.0 in the interior, smaller near edges.
fn init_damping(params: &SimParams, fields: &mut Fields) {
    let nx = params.nx;
    let nz = params.nz;
    let w = params.absorb_width as f32;
    let absorb_max = params.absorb_max;

    if params.absorb_width == 0 || absorb_max == 0.0 {
        fields.damping.fill(1.0);
        return;
    }

    for z in 0..nz {
        for x in 0..nx {
            let dist_x = (x as isize).min((nx - 1 - x) as isize) as f32;
            let dist_z = (z as isize).min((nz - 1 - z) as isize) as f32;
            let dist = dist_x.min(dist_z);

            let damp = if dist < w {
                let r = (w - dist) / w; // 0 at interior, 1 at boundary
                let sigma = absorb_max * r * r;
                1.0 - sigma
            } else {
                1.0
            };

            let id = idx(x, z, nx);
            fields.damping[id] = damp;
        }
    }
}

/// Ricker wavelet source.
fn ricker(t: f32, f0: f32) -> f32 {
    let pi = std::f32::consts::PI;
    let x = pi * f0 * (t - 3.0 / f0); // center pulse at t = 3/f0
    let x2 = x * x;
    (1.0 - 2.0 * x2) * (-x2).exp()
}

fn source_value(source: Source, step: usize, t: f32) -> f32 {
    match source {
        Source::Ricker { freq_hz } => ricker(t, freq_hz),
        Source::Impulse {
            step: source_step,
            amplitude,
        } => {
            if step == source_step {
                amplitude
            } else {
                0.0
            }
        }
    }
}

fn update_pressure(model: WaveModel, p: f32, p_prev: f32, lap: f32, coef_x: f32, dt: f32) -> f32 {
    match model {
        WaveModel::LosslessAcoustic => 2.0 * p - p_prev + coef_x * lap,
        WaveModel::LinearDampedAcoustic { gamma } => {
            let half_step_damping = 0.5 * gamma * dt;
            (2.0 * p - (1.0 - half_step_damping) * p_prev + coef_x * lap)
                / (1.0 + half_step_damping)
        }
    }
}

fn run_simulation(params: &SimParams, probes: &[Probe], log_progress: bool) -> SimOutput {
    assert_eq!(
        params.dx, params.dz,
        "the current 5-point Laplacian assumes dx == dz"
    );
    assert!(
        params.source_x > 0
            && params.source_x + 1 < params.nx
            && params.source_z > 0
            && params.source_z + 1 < params.nz,
        "source must be inside the finite-difference interior"
    );
    for probe in probes {
        assert!(probe.x < params.nx && probe.z < params.nz);
    }

    let nx = params.nx;
    let nz = params.nz;
    let n = nx * nz;

    let mut fields = Fields::new(nx, nz);
    init_damping(params, &mut fields);

    let c = params.c0;
    let dt = params.dt;
    let dx = params.dx;
    let coef = (c * dt / dx).powi(2);
    let source_idx = idx(params.source_x, params.source_z, nx);
    let mut probe_traces = vec![Vec::with_capacity(params.n_steps); probes.len()];

    for step in 0..params.n_steps {
        let t = step as f32 * dt;

        for v in fields.p_next.iter_mut() {
            *v = 0.0;
        }

        for z in 1..(nz - 1) {
            for x in 1..(nx - 1) {
                let i = idx(x, z, nx);
                let p = fields.p[i];

                let pxm = fields.p[idx(x - 1, z, nx)];
                let pxp = fields.p[idx(x + 1, z, nx)];
                let pzm = fields.p[idx(x, z - 1, nx)];
                let pzp = fields.p[idx(x, z + 1, nx)];

                let lap = pxm + pxp + pzm + pzp - 4.0 * p;

                fields.p_next[i] =
                    update_pressure(params.model, p, fields.p_prev[i], lap, coef, dt);
            }
        }

        fields.p_next[source_idx] += source_value(params.source, step, t);

        for i in 0..n {
            let damp = fields.damping[i];
            fields.p_next[i] *= damp;
            fields.p[i] *= damp;
            fields.p_prev[i] *= damp;
        }

        std::mem::swap(&mut fields.p_prev, &mut fields.p);
        std::mem::swap(&mut fields.p, &mut fields.p_next);

        for (trace, probe) in probe_traces.iter_mut().zip(probes.iter()) {
            trace.push(fields.p[idx(probe.x, probe.z, nx)]);
        }

        if log_progress && step % 100 == 0 {
            println!("Step {}/{}", step, params.n_steps);
        }
    }

    SimOutput {
        final_pressure: fields.p,
        probe_traces,
    }
}

/// Run the simulation and write a CSV of final pressure.
fn run_sim(params: &SimParams) {
    let output = run_simulation(params, &[], true);
    // Write final pressure field to CSV for plotting
    save_pressure_csv("output/pressure_final.csv", params, &output.final_pressure)
        .expect("Failed to write CSV");
    println!("Saved final pressure to output/pressure_final.csv");
}

fn save_pressure_csv(path: &str, params: &SimParams, p: &[f32]) -> std::io::Result<()> {
    create_dir_all("output")?;
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Each row is one z, comma-separated p(x,z)
    for z in 0..params.nz {
        for x in 0..params.nx {
            let i = idx(x, z, params.nx);
            write!(writer, "{:.6}", p[i])?;
            if x + 1 < params.nx {
                write!(writer, ",")?;
            }
        }
        writeln!(writer)?;
    }

    Ok(())
}

fn main() {
    let params = SimParams::default();
    run_sim(&params);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fast_baseline_params(model: WaveModel) -> SimParams {
        let nx = 61;
        let nz = 61;
        let dx = 0.05;
        let dz = 0.05;
        let c0 = 1.0;
        let dt = 0.015;

        SimParams {
            nx,
            nz,
            dx,
            dz,
            dt,
            n_steps: 150,
            c0,
            source_x: nx / 2,
            source_z: nz / 2,
            source: Source::Ricker { freq_hz: 5.0 },
            absorb_width: 0,
            absorb_max: 0.0,
            model,
        }
    }

    fn peak_time(trace: &[f32], dt: f32) -> (usize, f32, f32) {
        let (peak_step, peak_value) = trace
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.abs().total_cmp(&b.abs()))
            .map(|(step, value)| (step, *value))
            .expect("probe trace must not be empty");
        (peak_step, peak_step as f32 * dt, peak_value.abs())
    }

    fn trace_energy(trace: &[f32]) -> f32 {
        trace.iter().map(|value| value * value).sum()
    }

    #[test]
    fn impulse_source_only_fires_on_selected_step() {
        let source = Source::Impulse {
            step: 3,
            amplitude: 2.5,
        };

        assert_eq!(source_value(source, 2, 0.0), 0.0);
        assert_eq!(source_value(source, 3, 0.0), 2.5);
        assert_eq!(source_value(source, 4, 0.0), 0.0);
    }

    #[test]
    fn lossless_acoustic_pulse_arrives_near_expected_time() {
        let params = fast_baseline_params(WaveModel::LosslessAcoustic);
        let distance = 1.0;
        let probe = Probe {
            x: params.source_x + (distance / params.dx) as usize,
            z: params.source_z,
        };

        let output = run_simulation(&params, &[probe], false);
        let (_, observed_peak_time, peak_amp) = peak_time(&output.probe_traces[0], params.dt);
        let source_peak_time = 3.0 / 5.0;
        let expected_peak_time = source_peak_time + distance / params.c0;

        assert!(
            (observed_peak_time - expected_peak_time).abs() < 0.25,
            "lossless pulse peak arrived at {observed_peak_time}, expected near {expected_peak_time}"
        );
        assert!(peak_amp > 0.01, "probe should see a measurable pulse");
    }

    #[test]
    fn linear_damped_acoustic_keeps_arrival_time_but_loses_energy() {
        let lossless = fast_baseline_params(WaveModel::LosslessAcoustic);
        let damped = fast_baseline_params(WaveModel::LinearDampedAcoustic { gamma: 0.5 });
        let distance = 1.0;
        let probe = Probe {
            x: lossless.source_x + (distance / lossless.dx) as usize,
            z: lossless.source_z,
        };

        let lossless_output = run_simulation(&lossless, &[probe], false);
        let damped_output = run_simulation(&damped, &[probe], false);

        let (_, lossless_peak_time, lossless_peak_amp) =
            peak_time(&lossless_output.probe_traces[0], lossless.dt);
        let (_, damped_peak_time, damped_peak_amp) =
            peak_time(&damped_output.probe_traces[0], damped.dt);
        let lossless_energy = trace_energy(&lossless_output.probe_traces[0]);
        let damped_energy = trace_energy(&damped_output.probe_traces[0]);

        assert!(
            (damped_peak_time - lossless_peak_time).abs() < 0.08,
            "damped model should not materially change first-order travel time: lossless={lossless_peak_time}, damped={damped_peak_time}"
        );
        assert!(
            damped_peak_amp < lossless_peak_amp,
            "damped peak amplitude should be lower than lossless peak amplitude"
        );
        assert!(
            damped_energy < lossless_energy,
            "damped trace energy should be lower than lossless trace energy"
        );
    }
}
