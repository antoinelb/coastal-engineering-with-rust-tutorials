use std::error::Error;
use std::fmt;

#[derive(Debug)]
pub enum WaveError {
    InsufficientData,
}

impl fmt::Display for WaveError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            WaveError::InsufficientData => write!(f, "Insufficient wave data"),
        }
    }
}

impl Error for WaveError {}

#[derive(Debug, Clone)]
pub struct Wave {
    pub height: f64, // (m)
    pub period: f64, // (s)
}

#[derive(Debug)]
pub struct WaveMeasurement {
    pub time: f64,      // (s)
    pub elevation: f64, // (m)
}

pub struct WaveTimeSeries {
    pub measurements: Vec<WaveMeasurement>,
    pub sampling_rate: f64, // Hz
}

impl WaveTimeSeries {
    pub fn new(sampling_rate: f64) -> Self {
        WaveTimeSeries {
            measurements: Vec::new(),
            sampling_rate,
        }
    }

    pub fn add_measurement(&mut self, measurement: WaveMeasurement) {
        self.measurements.push(measurement);
    }

    pub fn mean_level(&self) -> Result<f64, WaveError> {
        if self.measurements.is_empty() {
            return Err(WaveError::InsufficientData);
        }

        let sum: f64 = self.measurements.iter().map(|m| m.elevation).sum();

        return Ok(sum / self.measurements.len() as f64);
    }

    pub fn detect_waves(&self) -> Result<Vec<Wave>, WaveError> {
        if self.measurements.len() < 3 {
            return Err(WaveError::InsufficientData);
        }

        let mean = self.mean_level()?;
        let detrended: Vec<f64> = self
            .measurements
            .iter()
            .map(|m| m.elevation - mean)
            .collect();

        let mut waves = Vec::new();
        let mut i = 0;

        while i < detrended.len() - 1 {
            // find next upward zero crossing
            if detrended[i] >= 0.0 && detrended[i + 1] < 0.0 {
                let start_idx = i;

                // find next downward zero crossig
                while i < detrended.len() - 1 && detrended[i + 1] < 0.0 {
                    i += 1;
                }

                if i < detrended.len() - 1 {
                    let end_idx = i;

                    let wave_segment = &self.measurements[start_idx..=end_idx];
                    // println!("{:?}", wave_segment);
                    let crest = wave_segment
                        .iter()
                        .map(|m| m.elevation)
                        .fold(f64::NEG_INFINITY, f64::max);
                    let trough = wave_segment
                        .iter()
                        .map(|m| m.elevation)
                        .fold(f64::INFINITY, f64::min);

                    let height = crest - trough;
                    let period = (end_idx - start_idx) as f64 / self.sampling_rate;

                    waves.push(Wave { height, period });
                }
            }

            i += 1;
        }

        if waves.is_empty() {
            return Err(WaveError::InsufficientData);
        }

        return Ok(waves);
    }

    pub fn calculate_statistics(&self) -> Result<WaveStatistics, WaveError> {
        let waves = self.detect_waves()?;

        if waves.is_empty() {
            return Err(WaveError::InsufficientData);
        }
        // println!("{:?}", waves);

        let mut sorted_waves = waves.clone();
        sorted_waves.sort_by(|a, b| b.height.partial_cmp(&a.height).unwrap());

        let n_third = (sorted_waves.len() as f64 / 3.0).ceil() as usize;
        let hs = sorted_waves[..n_third]
            .iter()
            .map(|w| w.height)
            .sum::<f64>()
            / n_third as f64;

        let hmax = sorted_waves[0].height;
        let wave_count = waves.len();
        let tmean = waves.iter().map(|w| w.period).sum::<f64>() / wave_count as f64;

        return Ok(WaveStatistics {
            hs,
            hmax,
            tmean,
            wave_count,
        });
    }

    pub fn assess_port_safety(&self, hs_threshold: f64) -> Result<bool, WaveError> {
        let stats = self.calculate_statistics()?;
        return Ok(stats.hs < hs_threshold);
    }

    pub fn rayleigh_parameters(&self) -> Result<(f64, f64), WaveError> {
        let waves = self.detect_waves()?;
        let hrms =
            (waves.iter().map(|w| w.height.powi(2)).sum::<f64>() / waves.len() as f64).sqrt();
        let scale = hrms / (std::f64::consts::PI / 2.0).sqrt();
        return Ok((hrms, scale));
    }

    pub fn exceedance_probability(&self, height: f64) -> Result<f64, WaveError> {
        let (_, scale) = self.rayleigh_parameters()?;
        return Ok((-height.powi(2) / (2.0 * scale.powi(2))).exp());
    }
}

#[derive(Debug)]
pub struct WaveStatistics {
    pub hs: f64,           // significant wave height
    pub hmax: f64,         // maximum wave height
    pub tmean: f64,        // mean wave period
    pub wave_count: usize, // number of waves
}

#[derive(Debug)]
pub struct SeasonalWaveClimate {
    pub winter_hs: f64,
    pub summer_hs: f64,
    pub storm_frequency: f64,
}
