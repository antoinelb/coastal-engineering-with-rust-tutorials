use crate::consts::G;
use crate::errors::WaveError;
use crate::waves::WaveTimeSeries;
use rustfft::{num_complex::Complex, FftPlanner};
use std::f64::consts::PI;

pub struct SpectralAnalysis {
    pub frequencies: Vec<f64>,
    pub spectrum: Vec<f64>,
    pub peak_frequency: f64,
    pub peak_period: f64,
}

pub struct WaveComponents {
    pub sea_spectrum: Vec<f64>,
    pub swell_spectrum: Vec<f64>,
    pub sea_hs: f64,
    pub swell_hs: f64,
    pub separation_frequency: f64,
}

pub struct DirectionalSpectrum {
    pub frequencies: Vec<f64>,
    pub directions: Vec<f64>,    // degrees
    pub spectrum: Vec<Vec<f64>>, // S(f, Θ)
}

impl WaveTimeSeries {
    // perform spectral analysis using fft
    pub fn perform_spectral_analysis(&self) -> Result<SpectralAnalysis, WaveError> {
        if self.measurements.len() < 256 {
            return Err(WaveError::InsufficientData);
        }

        // detrend data
        let mean = self.calculate_mean_level()?;
        let detrented: Vec<f64> = self
            .measurements
            .iter()
            .map(|m| m.elevation - mean)
            .collect();

        // apply hann window
        // reduces spectral leakage by tapering the signal to zero at both ends, improving
        // frequency
        let n = detrented.len();
        let windowed: Vec<Complex<f64>> = detrented
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                let window = 0.5 * (1.0 - (2.0 * PI * i as f64 / (n - 1) as f64).cos());
                return Complex::new(x * window, 0.0);
            })
            .collect();

        // perform fft
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(n);
        let mut spectrum_data = windowed.clone();
        fft.process(&mut spectrum_data);

        // calculate power spectral density
        let df = self.sampling_rate / n as f64;
        let frequencies: Vec<f64> = (0..n / 2).map(|i| i as f64 * df).collect();
        let spectrum: Vec<f64> = spectrum_data[..n / 2]
            .iter()
            .map(|c| 2.0 * c.norm_sqr() / (self.sampling_rate * n as f64))
            .collect();

        // find peak frequency
        let (peak_idx, _) = spectrum
            .iter()
            .enumerate()
            .skip(1) // skip DC component
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();
        let peak_frequency = frequencies[peak_idx];
        let peak_period = 1.0 / peak_frequency;

        return Ok(SpectralAnalysis {
            frequencies,
            spectrum,
            peak_frequency,
            peak_period,
        });
    }

    // calculate spectral moments
    pub fn spectral_moments(&self, order: i32) -> Result<f64, WaveError> {
        let analysis = self.perform_spectral_analysis()?;

        let moment = analysis
            .frequencies
            .iter()
            .zip(analysis.spectrum.iter())
            .skip(1) // skip DC component
            .map(|(f, s)| f.powi(order) * s * (analysis.frequencies[1] - analysis.frequencies[0]))
            .sum();

        return Ok(moment);
    }
}

impl SpectralAnalysis {
    // separate sea and swell components
    pub fn separate_sea_swell(&self, wind_speed: f64) -> WaveComponents {
        // Pierson-Moskowitz frequency for fully developed sea
        // f_m = 0.13 · g / (2 * pi * wind_speed)
        let fm = 0.13 * G / (2.0 * PI * wind_speed);

        // separation frequency (typical value)
        let separation_frequency = 1.2 * fm;

        // split spectrum
        let mut sea_spectrum = vec![0.0; self.spectrum.len()];
        let mut swell_spectrum = vec![0.0; self.spectrum.len()];
        for (i, (f, s)) in self
            .frequencies
            .iter()
            .zip(self.spectrum.iter())
            .enumerate()
        {
            if *f > separation_frequency {
                sea_spectrum[i] = *s;
            } else {
                swell_spectrum[i] = *s;
            }
        }

        // calculate wave heights
        let df = self.frequencies[1] - self.frequencies[0];
        let sea_m0 = sea_spectrum.iter().sum::<f64>() * df;
        let swell_m0 = swell_spectrum.iter().sum::<f64>() * df;
        let sea_hs = 4.0 * sea_m0.sqrt();
        let swell_hs = 4.0 * swell_m0.sqrt();

        return WaveComponents {
            sea_spectrum,
            swell_spectrum,
            sea_hs,
            swell_hs,
            separation_frequency,
        };
    }
}
