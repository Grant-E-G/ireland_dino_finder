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
