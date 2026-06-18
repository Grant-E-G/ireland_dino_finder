#![allow(dead_code)]

use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::path::Path;

use crate::grid::Grid;

#[derive(Clone, Debug)]
pub struct MaterialProperty {
    pub p_wave_velocity_m_s: Option<f32>,
}

#[derive(Clone, Debug)]
pub struct MaterialCatalog {
    by_id: HashMap<String, MaterialProperty>,
}

impl MaterialCatalog {
    pub fn from_properties_csv(path: impl AsRef<Path>) -> Result<Self, Box<dyn Error>> {
        let text = fs::read_to_string(path)?;
        let mut lines = text.lines();
        let header = lines.next().ok_or("material properties CSV is empty")?;
        let columns = split_csv_record(header);
        let id_col = find_column(&columns, "material_id")?;
        let velocity_col = find_column(&columns, "p_wave_velocity_m_s")?;

        let mut by_id = HashMap::new();
        for (line_no, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let fields = split_csv_record(line);
            if fields.len() <= id_col {
                return Err(format!("missing material_id on CSV line {}", line_no + 2).into());
            }
            let id = fields[id_col].trim().to_owned();
            let p_wave_velocity_m_s = fields
                .get(velocity_col)
                .and_then(|value| parse_optional_f32(value));

            by_id.insert(
                id.clone(),
                MaterialProperty {
                    p_wave_velocity_m_s,
                },
            );
        }

        Ok(Self { by_id })
    }

    pub fn velocity_for(&self, material_id: &str) -> Option<f32> {
        self.by_id
            .get(material_id)
            .and_then(|material| material.p_wave_velocity_m_s)
    }
}

#[derive(Clone, Debug)]
pub struct MaterialMap {
    pub grid: Grid,
    pub c: Vec<f32>,
    pub ids: Vec<String>,
}

impl MaterialMap {
    pub fn uniform(grid: Grid, material_id: impl Into<String>, c0: f32) -> Self {
        assert!(c0 > 0.0, "wave speed must be positive");
        let n = grid.len();
        let material_id = material_id.into();
        Self {
            grid,
            c: vec![c0; n],
            ids: vec![material_id; n],
        }
    }

    pub fn from_material_id_csv(
        path: impl AsRef<Path>,
        grid: Grid,
        catalog: &MaterialCatalog,
        fallback_velocity_m_s: f32,
    ) -> Result<Self, Box<dyn Error>> {
        assert!(
            fallback_velocity_m_s > 0.0,
            "fallback velocity must be positive"
        );
        let text = fs::read_to_string(path)?;
        let mut ids = Vec::with_capacity(grid.len());
        let mut c = Vec::with_capacity(grid.len());

        for line in text.lines().filter(|line| !line.trim().is_empty()) {
            for material_id in split_csv_record(line) {
                let trimmed = material_id.trim().to_owned();
                let velocity = catalog
                    .velocity_for(&trimmed)
                    .unwrap_or(fallback_velocity_m_s);
                if velocity <= 0.0 {
                    return Err(format!("material {trimmed} has non-positive wave speed").into());
                }
                ids.push(trimmed);
                c.push(velocity);
            }
        }

        if ids.len() != grid.len() {
            return Err(format!(
                "material map has {} cells, expected {} for {}x{} grid",
                ids.len(),
                grid.len(),
                grid.nx,
                grid.nz
            )
            .into());
        }

        Ok(Self { grid, c, ids })
    }

    pub fn max_speed(&self) -> f32 {
        self.c.iter().copied().fold(0.0_f32, f32::max)
    }
}

fn find_column(columns: &[String], name: &str) -> Result<usize, Box<dyn Error>> {
    columns
        .iter()
        .position(|column| column == name)
        .ok_or_else(|| format!("missing required CSV column {name}").into())
}

fn parse_optional_f32(value: &str) -> Option<f32> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        None
    } else {
        trimmed.parse().ok()
    }
}

pub fn split_csv_record(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut chars = line.chars().peekable();
    let mut in_quotes = false;

    while let Some(ch) = chars.next() {
        match ch {
            '"' if in_quotes && chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => {
                fields.push(field);
                field = String::new();
            }
            _ => field.push(ch),
        }
    }

    fields.push(field);
    fields
}
