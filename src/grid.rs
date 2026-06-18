pub fn idx(x: usize, z: usize, nx: usize) -> usize {
    z * nx + x
}

#[derive(Clone, Copy, Debug)]
pub struct Grid {
    pub nx: usize,
    pub nz: usize,
    pub dx: f32,
    pub dz: f32,
}

impl Grid {
    pub fn new(nx: usize, nz: usize, dx: f32, dz: f32) -> Self {
        assert!(nx >= 3, "grid nx must be at least 3");
        assert!(nz >= 3, "grid nz must be at least 3");
        assert!(dx > 0.0, "grid dx must be positive");
        assert!(dz > 0.0, "grid dz must be positive");
        Self { nx, nz, dx, dz }
    }

    pub fn len(&self) -> usize {
        self.nx * self.nz
    }

    pub fn id(&self, x: usize, z: usize) -> usize {
        idx(x, z, self.nx)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Probe {
    pub x: usize,
    pub z: usize,
}
