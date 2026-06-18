mod grid;
mod hdf5_helper;
mod materials;
mod model;
mod output;
mod solver;
mod source;

use solver::{OutputConfig, SimParams, run_simulation_with_output};

fn main() {
    let params = SimParams::default();
    run_simulation_with_output(&params, &[], true, Some(&OutputConfig::demo()))
        .expect("demo simulation failed");
}
