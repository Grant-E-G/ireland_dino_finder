mod grid;
mod hdf5_helper;
mod materials;
mod model;
mod output;
mod scenarios;
mod solver;
mod source;

use scenarios::{run_named_scenario, scenario_names};
use solver::{OutputConfig, SimParams, run_simulation_with_output};

fn main() {
    let mut args = std::env::args().skip(1);
    if let Some(flag_or_scenario) = args.next() {
        let scenario = if flag_or_scenario == "--scenario" {
            args.next().unwrap_or_else(|| {
                eprintln!(
                    "Missing scenario name. Available: {}",
                    scenario_names().join(", ")
                );
                std::process::exit(2);
            })
        } else {
            flag_or_scenario
        };

        if scenario == "--help" || scenario == "-h" {
            print_usage();
            return;
        }

        if let Err(err) = run_named_scenario(&scenario) {
            eprintln!("{err}");
            eprintln!("Available scenarios: {}", scenario_names().join(", "));
            std::process::exit(2);
        }
        return;
    }

    let params = SimParams::default();
    run_simulation_with_output(&params, &[], true, Some(&OutputConfig::demo()))
        .expect("demo simulation failed");
}

fn print_usage() {
    println!("Usage:");
    println!("  cargo run");
    println!("  cargo run -- --scenario <name>");
    println!();
    println!("Available scenarios:");
    for name in scenario_names() {
        println!("  {name}");
    }
}
