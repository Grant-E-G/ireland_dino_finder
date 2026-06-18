#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
pub enum Source {
    Ricker { freq_hz: f32 },
    Impulse { step: usize, amplitude: f32 },
}

/// Ricker wavelet source.
pub fn ricker(t: f32, f0: f32) -> f32 {
    let pi = std::f32::consts::PI;
    let x = pi * f0 * (t - 3.0 / f0);
    let x2 = x * x;
    (1.0 - 2.0 * x2) * (-x2).exp()
}

pub fn source_value(source: Source, step: usize, t: f32) -> f32 {
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
