use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::grid::Grid;

pub fn ensure_parent_dir(path: impl AsRef<Path>) -> std::io::Result<()> {
    if let Some(parent) = path.as_ref().parent() {
        if !parent.as_os_str().is_empty() {
            create_dir_all(parent)?;
        }
    }
    Ok(())
}

pub fn save_pressure_csv(path: impl AsRef<Path>, grid: Grid, p: &[f32]) -> std::io::Result<()> {
    ensure_parent_dir(&path)?;
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for z in 0..grid.nz {
        for x in 0..grid.nx {
            let i = grid.id(x, z);
            write!(writer, "{:.6}", p[i])?;
            if x + 1 < grid.nx {
                write!(writer, ",")?;
            }
        }
        writeln!(writer)?;
    }

    Ok(())
}

pub fn save_pressure_ppm(
    path: impl AsRef<Path>,
    grid: Grid,
    p: &[f32],
    scale: Option<f32>,
) -> std::io::Result<()> {
    ensure_parent_dir(&path)?;
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let vmax = scale.filter(|value| *value > 0.0).unwrap_or_else(|| {
        p.iter()
            .map(|value| value.abs())
            .fold(0.0_f32, f32::max)
            .max(1.0e-6)
    });

    writeln!(writer, "P3")?;
    writeln!(writer, "{} {}", grid.nx, grid.nz)?;
    writeln!(writer, "255")?;

    for z in (0..grid.nz).rev() {
        for x in 0..grid.nx {
            let value = (p[grid.id(x, z)] / vmax).clamp(-1.0, 1.0);
            let (r, g, b) = seismic_rgb(value);
            write!(writer, "{r} {g} {b} ")?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

fn seismic_rgb(value: f32) -> (u8, u8, u8) {
    if value >= 0.0 {
        let gb = ((1.0 - value) * 255.0).round() as u8;
        (255, gb, gb)
    } else {
        let rg = ((1.0 + value) * 255.0).round() as u8;
        (rg, rg, 255)
    }
}
