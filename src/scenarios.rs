use std::error::Error;

use crate::grid::Grid;
use crate::materials::MaterialMap;
use crate::model::WaveModel;
use crate::solver::{OutputConfig, SimParams, run_simulation_with_output};
use crate::source::Source;

pub fn run_named_scenario(name: &str) -> Result<(), Box<dyn Error>> {
    match name {
        "visual-all" => {
            for scenario in visual_scenarios() {
                run_visual_scenario(scenario)?;
            }
            Ok(())
        }
        _ => {
            let scenario = visual_scenarios()
                .into_iter()
                .find(|scenario| scenario.name == name)
                .ok_or_else(|| format!("unknown scenario '{name}'"))?;
            run_visual_scenario(scenario)
        }
    }
}

pub fn scenario_names() -> Vec<&'static str> {
    let mut names: Vec<&'static str> = visual_scenarios()
        .into_iter()
        .map(|scenario| scenario.name)
        .collect();
    names.insert(0, "visual-all");
    names
}

#[derive(Clone, Copy)]
struct VisualScenario {
    name: &'static str,
    model: WaveModel,
    material_case: MaterialCase,
}

#[derive(Clone, Copy)]
enum MaterialCase {
    Uniform,
    SlowBlock,
    FastTarget,
}

fn visual_scenarios() -> Vec<VisualScenario> {
    vec![
        VisualScenario {
            name: "visual_lossless_uniform",
            model: WaveModel::LosslessAcoustic,
            material_case: MaterialCase::Uniform,
        },
        VisualScenario {
            name: "visual_lossless_slow_block",
            model: WaveModel::LosslessAcoustic,
            material_case: MaterialCase::SlowBlock,
        },
        VisualScenario {
            name: "visual_lossless_fast_target",
            model: WaveModel::LosslessAcoustic,
            material_case: MaterialCase::FastTarget,
        },
        VisualScenario {
            name: "visual_damped_slow_block",
            model: WaveModel::LinearDampedAcoustic { gamma: 0.35 },
            material_case: MaterialCase::SlowBlock,
        },
        VisualScenario {
            name: "visual_sls_slow_block",
            model: WaveModel::StandardLinearSolid {
                damping_gamma: 0.35,
                relaxation_time_s: 0.10,
                relaxation_strength: 0.35,
            },
            material_case: MaterialCase::SlowBlock,
        },
        VisualScenario {
            name: "visual_fractional_constant_q_slow_block",
            model: WaveModel::FractionalConstantQ {
                q: 70.0,
                reference_freq_hz: 5.0,
                dispersion_strength: 0.20,
            },
            material_case: MaterialCase::SlowBlock,
        },
        VisualScenario {
            name: "visual_biot_sand_target",
            model: WaveModel::ReducedBiotPoroelastic {
                drag_gamma: 0.35,
                relaxation_time_s: 0.10,
                pore_coupling: 0.65,
            },
            material_case: MaterialCase::FastTarget,
        },
    ]
}

fn run_visual_scenario(scenario: VisualScenario) -> Result<(), Box<dyn Error>> {
    let params = visual_params(scenario);
    let output_dir = "output/visual";
    let config = OutputConfig {
        csv_path: Some(format!("{output_dir}/{}.csv", scenario.name)),
        hdf5_path: Some(format!("{output_dir}/{}.h5", scenario.name)),
        image_frame_dir: Some(format!("{output_dir}/{}_frames", scenario.name)),
        frame_interval: 2,
    };

    println!("Running scenario {}", scenario.name);
    run_simulation_with_output(&params, &[], true, Some(&config))?;
    println!(
        "Render with: python python/h5_to_mp4.py {output_dir}/{}.h5 {output_dir}/{}.mp4",
        scenario.name, scenario.name
    );
    println!(
        "No Python deps: ffmpeg -y -framerate 30 -i {output_dir}/{}_frames/frame_%06d.ppm -pix_fmt yuv420p {output_dir}/{}.mp4",
        scenario.name, scenario.name
    );
    Ok(())
}

fn visual_params(scenario: VisualScenario) -> SimParams {
    let nx = 120;
    let nz = 80;
    let dx = 0.025;
    let dz = 0.025;
    let c0 = 1.0;
    let grid = Grid::new(nx, nz, dx, dz);
    let material_map = visual_material_map(grid, scenario.material_case, c0);

    SimParams {
        grid,
        dt: 0.008,
        n_steps: 360,
        source_x: 15,
        source_z: nz / 2,
        source: Source::Ricker { freq_hz: 6.0 },
        absorb_width: 12,
        absorb_max: 0.04,
        model: scenario.model,
        material_map,
    }
}

fn visual_material_map(grid: Grid, material_case: MaterialCase, c0: f32) -> MaterialMap {
    let mut map = MaterialMap::uniform(grid, "background", c0);

    match material_case {
        MaterialCase::Uniform => {}
        MaterialCase::SlowBlock => {
            for z in 22..59 {
                for x in 52..77 {
                    let i = grid.id(x, z);
                    map.c[i] = 0.65;
                    map.ids[i] = "slow_block".to_owned();
                }
            }
        }
        MaterialCase::FastTarget => {
            for z in 16..65 {
                for x in 38..100 {
                    let i = grid.id(x, z);
                    map.c[i] = 0.78;
                    map.ids[i] = "sand_like_slow_background".to_owned();
                }
            }

            let cx = 70_i32;
            let cz = grid.nz as i32 / 2;
            let radius = 9_i32;
            for z in 0..grid.nz {
                for x in 0..grid.nx {
                    let dx = x as i32 - cx;
                    let dz = z as i32 - cz;
                    if dx * dx + dz * dz <= radius * radius {
                        let i = grid.id(x, z);
                        map.c[i] = 1.55;
                        map.ids[i] = "fast_bone_like_target".to_owned();
                    }
                }
            }
        }
    }

    map
}
