use hdf5_metno::dataset::Dataset;
use hdf5_metno::{Extents, SimpleExtents};
use hdf5_metno::{File, Result};
use ndarray::Array2;

pub struct Hdf5WaveWriter {
    file: File,
    ds_u: Dataset,
    nt_written: usize,
    ny: usize,
    nx: usize,
}

impl Hdf5WaveWriter {
    /// Create a new HDF5 file with an extendible dataset "u" of shape (time, ny, nx).
    ///
    /// - `path`: output .h5 file path
    /// - `ny`, `nx`: spatial grid size
    /// - `chunk_t`: chunk size along the time axis (e.g. 8, 16, 32)
    pub fn create(path: &str, ny: usize, nx: usize, chunk_t: usize) -> Result<Self> {
        // Create (or overwrite) the file
        let file = File::create(path)?;

        // Current shape = (0, ny, nx); all dims resizable (time will grow)
        let extents = SimpleExtents::resizable(vec![0, ny, nx]);
        let shape: Extents = extents.into();

        // NOTE: depending on crate version this may be `new_dataset::<T>()`
        // rather than `new_dataset_builder::<T>()`. Adjust if needed.
        let ds_u = file
            .new_dataset::<f32>() // dataset element type
            .shape(shape) // current + max dimensions
            .chunk((chunk_t, ny, nx)) // HDF5 chunk shape
            .create("u")?; // dataset name

        Ok(Self {
            file,
            ds_u,
            nt_written: 0,
            ny,
            nx,
        })
    }

    /// Optionally store coordinates and metadata once at the start.
    pub fn write_metadata(
        &self,
        x: Option<&[f64]>,
        y: Option<&[f64]>,
        dt: Option<f64>,
        dx: Option<f64>,
        dy: Option<f64>,
    ) -> Result<()> {
        if let Some(xvals) = x {
            self.file
                .new_dataset::<f64>()
                .shape(xvals.len())
                .create("x")?
                .write(xvals)?;
        }

        if let Some(yvals) = y {
            self.file
                .new_dataset::<f64>()
                .shape(yvals.len())
                .create("y")?
                .write(yvals)?;
        }

        if let Some(dt) = dt {
            self.file
                .new_attr::<f64>()
                .create("dt")?
                .write_scalar(&dt)?;
        }

        if let Some(dx) = dx {
            self.file
                .new_attr::<f64>()
                .create("dx")?
                .write_scalar(&dx)?;
        }

        if let Some(dy) = dy {
            self.file
                .new_attr::<f64>()
                .create("dy")?
                .write_scalar(&dy)?;
        }

        Ok(())
    }

    /// Append a single time slice u_t(y, x) to the "u" dataset.
    pub fn append_timestep(&mut self, u: &Array2<f32>) -> Result<()> {
        assert_eq!(u.shape(), &[self.ny, self.nx]);

        let cur_t = self.nt_written;
        let new_t = cur_t + 1;

        // Grow dataset from (cur_t, ny, nx) to (new_t, ny, nx)
        self.ds_u.resize((new_t, self.ny, self.nx))?;

        // Write directly as a 2D slice at fixed time index `cur_t`
        self.ds_u
            .write_slice(u.view(), (cur_t, 0..self.ny, 0..self.nx))?;

        self.nt_written = new_t;
        Ok(())
    }
}
