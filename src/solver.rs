use std::collections::HashMap;
use std::error::Error;

use ndarray::Array2;

use crate::grid::{Grid, Probe};
use crate::hdf5_helper::Hdf5WaveWriter;
use crate::materials::MaterialMap;
use crate::model::{WaveModel, update_pressure};
use crate::output::{ensure_parent_dir, save_pressure_csv, save_pressure_ppm};
use crate::source::{Source, source_value};

pub struct SimParams {
    pub grid: Grid,
    pub dt: f32,
    pub n_steps: usize,
    pub source_x: usize,
    pub source_z: usize,
    pub source: Source,
    pub absorb_width: usize,
    pub absorb_max: f32,
    pub model: WaveModel,
    pub material_map: MaterialMap,
}

impl SimParams {
    pub fn default() -> Self {
        let nx = 200;
        let nz = 200;
        let dx = 1.0;
        let dz = 1.0;
        let c0 = 2000.0;
        let dt = 0.00025;
        let grid = Grid::new(nx, nz, dx, dz);
        let material_map = MaterialMap::uniform(grid, "fresh_water_20c", c0);

        let r = c0 * dt / dx;
        println!("Courant number r = {}", r);

        Self {
            grid,
            dt,
            n_steps: 2000,
            source_x: nx / 2,
            source_z: nz / 2,
            source: Source::Ricker { freq_hz: 100.0 },
            absorb_width: 20,
            absorb_max: 0.015,
            model: WaveModel::LosslessAcoustic,
            material_map,
        }
    }
}

pub struct SimOutput {
    #[allow(dead_code)]
    pub final_pressure: Vec<f32>,
    #[allow(dead_code)]
    pub probe_traces: Vec<Vec<f32>>,
}

#[derive(Clone, Debug)]
pub struct OutputConfig {
    pub csv_path: Option<String>,
    pub hdf5_path: Option<String>,
    pub image_frame_dir: Option<String>,
    pub frame_interval: usize,
}

impl OutputConfig {
    pub fn demo() -> Self {
        Self {
            csv_path: Some("output/pressure_final.csv".to_owned()),
            hdf5_path: Some("output/wavefield.h5".to_owned()),
            image_frame_dir: None,
            frame_interval: 20,
        }
    }
}

struct Fields {
    p: Vec<f32>,
    p_prev: Vec<f32>,
    p_next: Vec<f32>,
    aux: Vec<f32>,
    aux_next: Vec<f32>,
    damping: Vec<f32>,
}

impl Fields {
    fn new(grid: Grid) -> Self {
        let n = grid.len();
        Self {
            p: vec![0.0; n],
            p_prev: vec![0.0; n],
            p_next: vec![0.0; n],
            aux: vec![0.0; n],
            aux_next: vec![0.0; n],
            damping: vec![1.0; n],
        }
    }
}

#[allow(dead_code)]
pub fn run_simulation(params: &SimParams, probes: &[Probe], log_progress: bool) -> SimOutput {
    run_simulation_with_output(params, probes, log_progress, None)
        .expect("simulation without output should not fail")
}

pub fn run_simulation_with_output(
    params: &SimParams,
    probes: &[Probe],
    log_progress: bool,
    output_config: Option<&OutputConfig>,
) -> Result<SimOutput, Box<dyn Error>> {
    validate_params(params, probes);

    let grid = params.grid;
    let n = grid.len();

    let mut fields = Fields::new(grid);
    init_damping(params, &mut fields);

    let dt = params.dt;
    let dx = grid.dx;
    let source_idx = grid.id(params.source_x, params.source_z);
    let mut probe_traces = vec![Vec::with_capacity(params.n_steps); probes.len()];
    let mut hdf5_writer = create_hdf5_writer(output_config, grid, params)?;

    for step in 0..params.n_steps {
        let t = step as f32 * dt;

        for v in fields.p_next.iter_mut() {
            *v = 0.0;
        }
        fields.aux_next.copy_from_slice(&fields.aux);

        for z in 1..(grid.nz - 1) {
            for x in 1..(grid.nx - 1) {
                let i = grid.id(x, z);
                let p = fields.p[i];

                let pxm = fields.p[grid.id(x - 1, z)];
                let pxp = fields.p[grid.id(x + 1, z)];
                let pzm = fields.p[grid.id(x, z - 1)];
                let pzp = fields.p[grid.id(x, z + 1)];

                let lap = pxm + pxp + pzm + pzp - 4.0 * p;
                let c = params.material_map.c[i];
                let coef = (c * dt / dx).powi(2);

                let (p_next, aux_next) = update_pressure(
                    params.model,
                    p,
                    fields.p_prev[i],
                    lap,
                    fields.aux[i],
                    coef,
                    dt,
                );
                fields.p_next[i] = p_next;
                fields.aux_next[i] = aux_next;
            }
        }

        fields.p_next[source_idx] += source_value(params.source, step, t);

        for i in 0..n {
            let damp = fields.damping[i];
            fields.p_next[i] *= damp;
            fields.p[i] *= damp;
            fields.p_prev[i] *= damp;
            fields.aux_next[i] *= damp;
            fields.aux[i] *= damp;
        }

        std::mem::swap(&mut fields.p_prev, &mut fields.p);
        std::mem::swap(&mut fields.p, &mut fields.p_next);
        std::mem::swap(&mut fields.aux, &mut fields.aux_next);

        for (trace, probe) in probe_traces.iter_mut().zip(probes.iter()) {
            trace.push(fields.p[grid.id(probe.x, probe.z)]);
        }

        write_frame_if_needed(&mut hdf5_writer, output_config, step, grid, &fields.p)?;
        write_image_frame_if_needed(output_config, step, grid, &fields.p)?;

        if log_progress && step % 100 == 0 {
            println!("Step {}/{}", step, params.n_steps);
        }
    }

    if let Some(config) = output_config {
        if let Some(path) = &config.csv_path {
            save_pressure_csv(path, grid, &fields.p)?;
            println!("Saved final pressure to {path}");
        }
    }

    write_probe_metadata_if_needed(&hdf5_writer, probes, &probe_traces)?;

    Ok(SimOutput {
        final_pressure: fields.p,
        probe_traces,
    })
}

fn validate_params(params: &SimParams, probes: &[Probe]) {
    assert_eq!(
        params.grid.dx, params.grid.dz,
        "the current 5-point Laplacian assumes dx == dz"
    );
    assert_eq!(
        params.material_map.grid.len(),
        params.grid.len(),
        "material map size must match simulation grid"
    );
    let courant = params.material_map.max_speed() * params.dt / params.grid.dx;
    assert!(
        courant < 1.0 / 2.0_f32.sqrt(),
        "unstable Courant number {courant}; expected < 1/sqrt(2)"
    );
    assert!(
        params.source_x > 0
            && params.source_x + 1 < params.grid.nx
            && params.source_z > 0
            && params.source_z + 1 < params.grid.nz,
        "source must be inside the finite-difference interior"
    );
    for probe in probes {
        assert!(probe.x < params.grid.nx && probe.z < params.grid.nz);
    }
}

fn init_damping(params: &SimParams, fields: &mut Fields) {
    let grid = params.grid;
    let w = params.absorb_width as f32;
    let absorb_max = params.absorb_max;

    if params.absorb_width == 0 || absorb_max == 0.0 {
        fields.damping.fill(1.0);
        return;
    }

    for z in 0..grid.nz {
        for x in 0..grid.nx {
            let dist_x = (x as isize).min((grid.nx - 1 - x) as isize) as f32;
            let dist_z = (z as isize).min((grid.nz - 1 - z) as isize) as f32;
            let dist = dist_x.min(dist_z);

            let damp = if dist < w {
                let r = (w - dist) / w;
                let sigma = absorb_max * r * r;
                1.0 - sigma
            } else {
                1.0
            };

            let id = grid.id(x, z);
            fields.damping[id] = damp;
        }
    }
}

fn create_hdf5_writer(
    output_config: Option<&OutputConfig>,
    grid: Grid,
    params: &SimParams,
) -> Result<Option<Hdf5WaveWriter>, Box<dyn Error>> {
    let Some(config) = output_config else {
        return Ok(None);
    };
    let Some(path) = &config.hdf5_path else {
        return Ok(None);
    };
    if config.frame_interval == 0 {
        return Err("HDF5 frame_interval must be greater than zero".into());
    }

    ensure_parent_dir(path)?;
    let writer = Hdf5WaveWriter::create(path, grid.nz, grid.nx, 8)?;
    let x: Vec<f64> = (0..grid.nx).map(|x| x as f64 * grid.dx as f64).collect();
    let y: Vec<f64> = (0..grid.nz).map(|z| z as f64 * grid.dz as f64).collect();
    writer.write_metadata(
        Some(&x),
        Some(&y),
        Some(params.dt as f64),
        Some(grid.dx as f64),
        Some(grid.dz as f64),
    )?;
    writer.write_static_field_f32("material_speed", &params.material_map.c)?;
    let (material_indices, material_name_bytes, material_count, material_name_width) =
        material_id_metadata(&params.material_map.ids);
    writer.write_static_field_i32("material_id_index", &material_indices)?;
    writer.write_byte_table(
        "material_id_names",
        &material_name_bytes,
        material_count,
        material_name_width,
    )?;
    let source_wavelet: Vec<f32> = (0..params.n_steps)
        .map(|step| source_value(params.source, step, step as f32 * params.dt))
        .collect();
    writer.write_vector_f32("source_wavelet", &source_wavelet)?;
    Ok(Some(writer))
}

fn material_id_metadata(ids: &[String]) -> (Vec<i32>, Vec<u8>, usize, usize) {
    let mut by_id: HashMap<&str, i32> = HashMap::new();
    let mut names: Vec<&str> = Vec::new();
    let mut indices = Vec::with_capacity(ids.len());

    for id in ids {
        let next_index = names.len() as i32;
        let index = *by_id.entry(id.as_str()).or_insert_with(|| {
            names.push(id.as_str());
            next_index
        });
        indices.push(index);
    }

    let name_width = names
        .iter()
        .map(|name| name.len())
        .max()
        .unwrap_or(1)
        .max(1);
    let mut name_bytes = vec![0_u8; names.len().max(1) * name_width];
    for (row, name) in names.iter().enumerate() {
        let offset = row * name_width;
        let bytes = name.as_bytes();
        name_bytes[offset..offset + bytes.len()].copy_from_slice(bytes);
    }

    (indices, name_bytes, names.len().max(1), name_width)
}

fn write_probe_metadata_if_needed(
    hdf5_writer: &Option<Hdf5WaveWriter>,
    probes: &[Probe],
    probe_traces: &[Vec<f32>],
) -> Result<(), Box<dyn Error>> {
    let Some(writer) = hdf5_writer else {
        return Ok(());
    };
    if probes.is_empty() {
        return Ok(());
    }

    let probe_x: Vec<i32> = probes.iter().map(|probe| probe.x as i32).collect();
    let probe_z: Vec<i32> = probes.iter().map(|probe| probe.z as i32).collect();
    writer.write_vector_i32("probe_x", &probe_x)?;
    writer.write_vector_i32("probe_z", &probe_z)?;

    let n_probes = probe_traces.len();
    let n_steps = probe_traces.first().map_or(0, Vec::len);
    let mut flattened = Vec::with_capacity(n_probes * n_steps);
    for trace in probe_traces {
        assert_eq!(
            trace.len(),
            n_steps,
            "all probe traces must have same length"
        );
        flattened.extend_from_slice(trace);
    }
    writer.write_matrix_f32("probe_traces", &flattened, n_probes, n_steps)?;
    Ok(())
}

fn write_frame_if_needed(
    hdf5_writer: &mut Option<Hdf5WaveWriter>,
    output_config: Option<&OutputConfig>,
    step: usize,
    grid: Grid,
    pressure: &[f32],
) -> Result<(), Box<dyn Error>> {
    let Some(writer) = hdf5_writer else {
        return Ok(());
    };
    let Some(config) = output_config else {
        return Ok(());
    };
    if step % config.frame_interval != 0 {
        return Ok(());
    }

    let frame = Array2::from_shape_vec((grid.nz, grid.nx), pressure.to_vec())?;
    writer.append_timestep(&frame)?;
    Ok(())
}

fn write_image_frame_if_needed(
    output_config: Option<&OutputConfig>,
    step: usize,
    grid: Grid,
    pressure: &[f32],
) -> Result<(), Box<dyn Error>> {
    let Some(config) = output_config else {
        return Ok(());
    };
    let Some(dir) = &config.image_frame_dir else {
        return Ok(());
    };
    if config.frame_interval == 0 {
        return Err("image frame_interval must be greater than zero".into());
    }
    if step % config.frame_interval != 0 {
        return Ok(());
    }

    std::fs::create_dir_all(dir)?;
    let frame_index = step / config.frame_interval;
    let path = format!("{dir}/frame_{frame_index:06}.ppm");
    save_pressure_ppm(path, grid, pressure, None)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;
    use crate::materials::{MaterialCatalog, MaterialMap, split_csv_record};
    use crate::model::{WaveModel, grunwald_letnikov_weights};
    use crate::source::{Source, ricker, source_value};

    fn fast_baseline_params(model: WaveModel) -> SimParams {
        let nx = 61;
        let nz = 61;
        let dx = 0.05;
        let dz = 0.05;
        let c0 = 1.0;
        let dt = 0.015;
        let grid = Grid::new(nx, nz, dx, dz);

        SimParams {
            grid,
            dt,
            n_steps: 150,
            source_x: nx / 2,
            source_z: nz / 2,
            source: Source::Ricker { freq_hz: 5.0 },
            absorb_width: 0,
            absorb_max: 0.0,
            model,
            material_map: MaterialMap::uniform(grid, "test_medium", c0),
        }
    }

    fn fast_baseline_probe(params: &SimParams) -> Probe {
        let distance = 1.0;
        Probe {
            x: params.source_x + (distance / params.grid.dx) as usize,
            z: params.source_z,
        }
    }

    fn finite_values(values: &[f32]) -> bool {
        values.iter().all(|value| value.is_finite())
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

    fn field_energy(field: &[f32]) -> f32 {
        field.iter().map(|value| value * value).sum()
    }

    fn assert_lossy_model_is_reasonable(model: WaveModel) {
        let lossless = fast_baseline_params(WaveModel::LosslessAcoustic);
        let candidate = fast_baseline_params(model);
        let probe = fast_baseline_probe(&lossless);

        let lossless_output = run_simulation(&lossless, &[probe], false);
        let candidate_output = run_simulation(&candidate, &[probe], false);

        let (_, lossless_peak_time, lossless_peak_amp) =
            peak_time(&lossless_output.probe_traces[0], lossless.dt);
        let (_, candidate_peak_time, candidate_peak_amp) =
            peak_time(&candidate_output.probe_traces[0], candidate.dt);
        let lossless_energy = trace_energy(&lossless_output.probe_traces[0]);
        let candidate_energy = trace_energy(&candidate_output.probe_traces[0]);

        assert!(
            finite_values(&candidate_output.final_pressure),
            "model produced non-finite final pressure values: {model:?}"
        );
        assert!(
            finite_values(&candidate_output.probe_traces[0]),
            "model produced non-finite probe trace values: {model:?}"
        );
        assert!(
            (candidate_peak_time - lossless_peak_time).abs() < 0.12,
            "lossy model should preserve first-order travel time: model={model:?}, lossless={lossless_peak_time}, candidate={candidate_peak_time}"
        );
        assert!(
            candidate_peak_amp < lossless_peak_amp,
            "lossy model should lower peak amplitude: model={model:?}, lossless={lossless_peak_amp}, candidate={candidate_peak_amp}"
        );
        assert!(
            candidate_energy < lossless_energy,
            "lossy model should lower probe energy: model={model:?}, lossless={lossless_energy}, candidate={candidate_energy}"
        );
    }

    fn temp_path(name: &str) -> String {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "ireland_dino_finder_{}_{}",
            std::process::id(),
            name
        ));
        path.to_string_lossy().into_owned()
    }

    #[test]
    fn csv_splitter_handles_quoted_commas() {
        let fields = split_csv_record("a,\"b, c\",\"d\"\"e\"");
        assert_eq!(fields, vec!["a", "b, c", "d\"e"]);
    }

    #[test]
    fn material_catalog_and_id_map_load_from_csv() {
        let properties_path = temp_path("materials.csv");
        let map_path = temp_path("map.csv");
        fs::write(
            &properties_path,
            "material_id,p_wave_velocity_m_s,notes\nwater,1.0,\"quoted, note\"\nbone,2.0,\n",
        )
        .unwrap();
        fs::write(
            &map_path,
            "water,bone,water\nbone,water,bone\nwater,bone,water\n",
        )
        .unwrap();

        let catalog = MaterialCatalog::from_properties_csv(&properties_path).unwrap();
        let grid = Grid::new(3, 3, 1.0, 1.0);
        let map = MaterialMap::from_material_id_csv(&map_path, grid, &catalog, 1.5).unwrap();

        assert_eq!(map.c, vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0]);
        assert_eq!(map.ids[1], "bone");

        let _ = fs::remove_file(properties_path);
        let _ = fs::remove_file(map_path);
    }

    #[test]
    fn ppm_mask_loads_material_ids_from_colors() {
        let properties_path = temp_path("ppm_materials.csv");
        let mask_path = temp_path("mask.ppm");
        fs::write(
            &properties_path,
            "material_id,p_wave_velocity_m_s\nwater,1.0\nbone,2.0\n",
        )
        .unwrap();
        fs::write(
            &mask_path,
            "P3\n# water black, bone white\n3 3\n255\n0 0 0 255 255 255 0 0 0\n255 255 255 0 0 0 255 255 255\n0 0 0 255 255 255 0 0 0\n",
        )
        .unwrap();

        let catalog = MaterialCatalog::from_properties_csv(&properties_path).unwrap();
        let grid = Grid::new(3, 3, 1.0, 1.0);
        let map = MaterialMap::from_ppm_mask(
            &mask_path,
            grid,
            &[((0, 0, 0), "water"), ((255, 255, 255), "bone")],
            &catalog,
            1.5,
        )
        .unwrap();

        assert_eq!(
            map.ids,
            vec![
                "water", "bone", "water", "bone", "water", "bone", "water", "bone", "water"
            ]
        );
        assert_eq!(map.c, vec![1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0]);

        let _ = fs::remove_file(properties_path);
        let _ = fs::remove_file(mask_path);
    }

    #[test]
    fn ricker_source_peaks_near_center_time() {
        let f0 = 5.0;
        let center = 3.0 / f0;
        let peak = ricker(center, f0);
        let before = ricker(center - 0.03, f0);
        let after = ricker(center + 0.03, f0);

        assert!((peak - 1.0).abs() < 1.0e-6);
        assert!(peak > before);
        assert!(peak > after);
    }

    #[test]
    fn grunwald_weights_start_with_expected_terms() {
        let weights = grunwald_letnikov_weights(0.5, 4);

        assert!((weights[0] - 1.0).abs() < 1.0e-6);
        assert!((weights[1] + 0.5).abs() < 1.0e-6);
        assert!((weights[2] + 0.125).abs() < 1.0e-6);
        assert!((weights[3] + 0.0625).abs() < 1.0e-6);
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
        let probe = fast_baseline_probe(&params);

        let output = run_simulation(&params, &[probe], false);
        let (_, observed_peak_time, peak_amp) = peak_time(&output.probe_traces[0], params.dt);
        let source_peak_time = 3.0 / 5.0;
        let expected_peak_time = source_peak_time + distance / 1.0;

        assert!(
            (observed_peak_time - expected_peak_time).abs() < 0.25,
            "lossless pulse peak arrived at {observed_peak_time}, expected near {expected_peak_time}"
        );
        assert!(peak_amp > 0.01, "probe should see a measurable pulse");
        assert!(finite_values(&output.final_pressure));
        assert!(finite_values(&output.probe_traces[0]));
    }

    #[test]
    fn heterogeneous_speed_map_changes_travel_time() {
        let uniform = fast_baseline_params(WaveModel::LosslessAcoustic);
        let mut heterogeneous = fast_baseline_params(WaveModel::LosslessAcoustic);
        for z in 0..heterogeneous.grid.nz {
            for x in (heterogeneous.source_x + 1)..heterogeneous.grid.nx {
                let i = heterogeneous.grid.id(x, z);
                heterogeneous.material_map.c[i] = 0.7;
                heterogeneous.material_map.ids[i] = "slow_zone".to_owned();
            }
        }
        let probe = fast_baseline_probe(&uniform);

        let uniform_output = run_simulation(&uniform, &[probe], false);
        let heterogeneous_output = run_simulation(&heterogeneous, &[probe], false);
        let (_, uniform_peak_time, _) = peak_time(&uniform_output.probe_traces[0], uniform.dt);
        let (_, heterogeneous_peak_time, _) =
            peak_time(&heterogeneous_output.probe_traces[0], heterogeneous.dt);

        assert!(
            heterogeneous_peak_time > uniform_peak_time,
            "slower material path should delay the probe peak"
        );
    }

    #[test]
    fn linear_damped_acoustic_keeps_arrival_time_but_loses_energy() {
        let lossless = fast_baseline_params(WaveModel::LosslessAcoustic);
        let damped = fast_baseline_params(WaveModel::LinearDampedAcoustic { gamma: 0.5 });
        let probe = fast_baseline_probe(&lossless);

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

    #[test]
    fn standard_linear_solid_baseline_is_reasonable() {
        assert_lossy_model_is_reasonable(WaveModel::StandardLinearSolid {
            damping_gamma: 0.5,
            relaxation_time_s: 0.12,
            relaxation_strength: 0.05,
        });
    }

    #[test]
    fn fractional_constant_q_baseline_is_reasonable() {
        assert_lossy_model_is_reasonable(WaveModel::FractionalConstantQ {
            q: 60.0,
            reference_freq_hz: 5.0,
            dispersion_strength: 0.05,
        });
    }

    #[test]
    fn reduced_biot_baseline_is_reasonable() {
        assert_lossy_model_is_reasonable(WaveModel::ReducedBiotPoroelastic {
            drag_gamma: 0.4,
            relaxation_time_s: 0.10,
            pore_coupling: 0.65,
        });
    }

    #[test]
    fn homogeneous_lossless_field_is_symmetric_about_source_axis() {
        let params = fast_baseline_params(WaveModel::LosslessAcoustic);
        let offset = 14;
        let left = Probe {
            x: params.source_x - offset,
            z: params.source_z,
        };
        let right = Probe {
            x: params.source_x + offset,
            z: params.source_z,
        };

        let output = run_simulation(&params, &[left, right], false);
        let max_delta = output.probe_traces[0]
            .iter()
            .zip(output.probe_traces[1].iter())
            .map(|(left_value, right_value)| (left_value - right_value).abs())
            .fold(0.0_f32, f32::max);

        assert!(
            max_delta < 1.0e-5,
            "symmetric probes diverged by {max_delta}"
        );
    }

    #[test]
    fn zero_source_keeps_field_at_rest() {
        let mut params = fast_baseline_params(WaveModel::LosslessAcoustic);
        params.source = Source::Impulse {
            step: params.n_steps + 1,
            amplitude: 1.0,
        };
        let probe = fast_baseline_probe(&params);

        let output = run_simulation(&params, &[probe], false);

        assert_eq!(field_energy(&output.final_pressure), 0.0);
        assert_eq!(trace_energy(&output.probe_traces[0]), 0.0);
    }

    #[test]
    fn sponge_boundary_reduces_late_field_energy() {
        let mut no_sponge = fast_baseline_params(WaveModel::LosslessAcoustic);
        no_sponge.n_steps = 260;

        let mut sponge = fast_baseline_params(WaveModel::LosslessAcoustic);
        sponge.n_steps = no_sponge.n_steps;
        sponge.absorb_width = 10;
        sponge.absorb_max = 0.08;

        let no_sponge_output = run_simulation(&no_sponge, &[], false);
        let sponge_output = run_simulation(&sponge, &[], false);

        assert!(
            field_energy(&sponge_output.final_pressure)
                < field_energy(&no_sponge_output.final_pressure),
            "sponge boundaries should reduce late field energy"
        );
    }

    #[test]
    fn hdf5_and_csv_output_are_written_from_active_solver_path() {
        let mut params = fast_baseline_params(WaveModel::LosslessAcoustic);
        params.n_steps = 8;
        let probe = fast_baseline_probe(&params);
        let csv_path = temp_path("pressure.csv");
        let hdf5_path = temp_path("wavefield.h5");
        let config = OutputConfig {
            csv_path: Some(csv_path.clone()),
            hdf5_path: Some(hdf5_path.clone()),
            image_frame_dir: None,
            frame_interval: 2,
        };

        run_simulation_with_output(&params, &[probe], false, Some(&config)).unwrap();

        assert!(fs::metadata(&csv_path).unwrap().len() > 0);
        assert!(fs::metadata(&hdf5_path).unwrap().len() > 0);

        let _ = fs::remove_file(csv_path);
        let _ = fs::remove_file(hdf5_path);
    }

    #[test]
    #[should_panic(expected = "unstable Courant number")]
    fn unstable_courant_number_is_rejected() {
        let mut params = fast_baseline_params(WaveModel::LosslessAcoustic);
        params.dt = 0.05;

        let _ = run_simulation(&params, &[], false);
    }
}
