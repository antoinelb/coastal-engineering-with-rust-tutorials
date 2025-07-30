use crate::errors::WaveError;

// a single wave measurement
#[derive(Debug, Clone)]
pub struct WaveMeasurement {
    pub time: f64,      // seconds since start
    pub elevation: f64, // water level in meters
}

// container for wave time series data
pub struct WaveTimeSeries {
    pub measurements: Vec<WaveMeasurement>,
    pub sampling_rate: f64, // in Hz
}

impl WaveTimeSeries {
    // creates a new wave time series
    pub fn new(sampling_rate: f64) -> Self {
        return WaveTimeSeries {
            measurements: Vec::new(),
            sampling_rate,
        };
    }

    // adds a new measurement
    pub fn add_measurement(&mut self, measurement: WaveMeasurement) {
        self.measurements.push(measurement);
    }

    // computes the mean water level
    pub fn calculate_mean_level(&self) -> Result<f64, WaveError> {
        if self.measurements.is_empty() {
            return Err(WaveError::InsufficientData);
        }

        let sum: f64 = self.measurements.iter().map(|m| m.elevation).sum();
        let mean: f64 = sum / self.measurements.len() as f64;

        return Ok(mean);
    }
}
