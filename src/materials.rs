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

    pub fn from_ppm_mask(
        path: impl AsRef<Path>,
        grid: Grid,
        color_materials: &[((u8, u8, u8), &str)],
        catalog: &MaterialCatalog,
        fallback_velocity_m_s: f32,
    ) -> Result<Self, Box<dyn Error>> {
        assert!(
            fallback_velocity_m_s > 0.0,
            "fallback velocity must be positive"
        );
        let text = fs::read_to_string(path)?;
        let mut tokens = ppm_tokens(&text);

        let magic = tokens.next().ok_or("PPM mask is empty")?;
        if magic != "P3" {
            return Err("only ASCII P3 PPM masks are supported".into());
        }

        let width = parse_ppm_usize(tokens.next(), "width")?;
        let height = parse_ppm_usize(tokens.next(), "height")?;
        let max_value = parse_ppm_usize(tokens.next(), "max value")?;
        if width != grid.nx || height != grid.nz {
            return Err(format!(
                "PPM mask is {}x{}, expected {}x{}",
                width, height, grid.nx, grid.nz
            )
            .into());
        }
        if max_value == 0 || max_value > 255 {
            return Err("PPM mask max value must be in 1..=255".into());
        }

        let mut by_color = HashMap::new();
        for (rgb, material_id) in color_materials {
            by_color.insert(*rgb, *material_id);
        }

        let mut ids = Vec::with_capacity(grid.len());
        let mut c = Vec::with_capacity(grid.len());
        for pixel_index in 0..grid.len() {
            let r = parse_ppm_color(tokens.next(), max_value, "red")?;
            let g = parse_ppm_color(tokens.next(), max_value, "green")?;
            let b = parse_ppm_color(tokens.next(), max_value, "blue")?;
            let material_id = by_color.get(&(r, g, b)).ok_or_else(|| {
                format!("unmapped PPM color #{r:02x}{g:02x}{b:02x} at pixel {pixel_index}")
            })?;
            let velocity = catalog
                .velocity_for(material_id)
                .unwrap_or(fallback_velocity_m_s);
            if velocity <= 0.0 {
                return Err(format!("material {material_id} has non-positive wave speed").into());
            }
            ids.push((*material_id).to_owned());
            c.push(velocity);
        }

        if tokens.next().is_some() {
            return Err("PPM mask has extra pixel data".into());
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

fn ppm_tokens(text: &str) -> impl Iterator<Item = &str> {
    text.lines()
        .map(|line| line.split_once('#').map_or(line, |(before, _)| before))
        .flat_map(str::split_whitespace)
}

fn parse_ppm_usize(value: Option<&str>, label: &str) -> Result<usize, Box<dyn Error>> {
    value
        .ok_or_else(|| format!("missing PPM {label}"))?
        .parse()
        .map_err(|_| format!("invalid PPM {label}").into())
}

fn parse_ppm_color(
    value: Option<&str>,
    max_value: usize,
    channel: &str,
) -> Result<u8, Box<dyn Error>> {
    let raw = parse_ppm_usize(value, channel)?;
    if raw > max_value {
        return Err(format!("PPM {channel} channel {raw} exceeds max value {max_value}").into());
    }
    Ok(((raw * 255 + max_value / 2) / max_value) as u8)
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
