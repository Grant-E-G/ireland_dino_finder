use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Write};

fn idx(x: usize, z: usize, nx: usize) -> usize {
    z * nx + x
}

struct SimParams {
    nx: usize,
    nz: usize,
    dx: f32,
    dz: f32,
    dt: f32,
    n_steps: usize,
    c0: f32,
    source_freq_hz: f32,
    source_duration_s: f32,
    absorb_width: usize,
    absorb_max: f32,
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
            source_freq_hz: 20.0, // Hz (very low; just for demo)
            source_duration_s: 0.05,
            absorb_width: 20,
            absorb_max: 0.015,
        }
    }
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

/// Run the simulation and write a CSV of final pressure.
fn run_sim(params: &SimParams) {
    let nx = params.nx;
    let nz = params.nz;
    let n = nx * nz;

    let mut fields = Fields::new(nx, nz);
    init_damping(params, &mut fields);

    // Precompute coefficient for Laplacian update (assuming dx = dz)
    let c = params.c0;
    let dt = params.dt;
    let dx = params.dx;
    let coef = (c * dt / dx).powi(2);

    // Source location: center of grid
    let sx = nx / 2;
    let sz = nz / 2;
    let source_idx = idx(sx, sz, nx);

    // Time-stepping loop
    for step in 0..params.n_steps {
        let t = step as f32 * dt;

        // Zero out next field
        for v in fields.p_next.iter_mut() {
            *v = 0.0;
        }

        // Interior update (2D 5-point Laplacian)
        for z in 1..(nz - 1) {
            for x in 1..(nx - 1) {
                let i = idx(x, z, nx);
                let p = fields.p[i];

                let pxm = fields.p[idx(x - 1, z, nx)];
                let pxp = fields.p[idx(x + 1, z, nx)];
                let pzm = fields.p[idx(x, z - 1, nx)];
                let pzp = fields.p[idx(x, z + 1, nx)];

                let lap = pxm + pxp + pzm + pzp - 4.0 * p;

                fields.p_next[i] = 2.0 * p - fields.p_prev[i] + coef * lap;
            }
        }

        // Add source term at center (in p_next)
        if t < params.source_duration_s {
            let s = ricker(t, params.source_freq_hz);
            fields.p_next[source_idx] += s;
        }

        // Apply damping (simple sponge)
        for i in 0..n {
            let damp = fields.damping[i];
            fields.p_next[i] *= damp;
            fields.p[i] *= damp;
            fields.p_prev[i] *= damp;
        }

        // Rotate time levels: p_prev <- p, p <- p_next
        std::mem::swap(&mut fields.p_prev, &mut fields.p);
        std::mem::swap(&mut fields.p, &mut fields.p_next);

        if step % 100 == 0 {
            println!("Step {}/{}", step, params.n_steps);
        }
    }

    // Write final pressure field to CSV for plotting
    save_pressure_csv("output/pressure_final.csv", params, &fields.p).expect("Failed to write CSV");
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
