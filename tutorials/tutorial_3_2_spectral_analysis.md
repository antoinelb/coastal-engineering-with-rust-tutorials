# Tutorial 3.2: Spectral Analysis Implementation
## Ocean Waves Chapter 3.5-3.6

### Learning Objectives
By the end of this tutorial, you will:
1. Understand wave energy spectra and their importance
2. Implement FFT-based spectral analysis in Rust
3. Calculate key spectral parameters (peak period, spectral moments)
4. Distinguish between sea and swell components
5. Apply spectral analysis to vessel motion prediction

### Prerequisites
- Completed Tutorial 3.1
- Read Coastal Dynamics Ch. 3.5-3.6
- Understanding of Fourier transforms

### Introduction: Why Spectral Analysis?

Wave fields are rarely composed of single periodic waves. Instead, they contain energy across multiple frequencies. Spectral analysis reveals:
- **Energy distribution**: Which periods contain most energy
- **Sea vs swell**: Locally generated vs distant waves
- **Directional spreading**: Wave energy by direction
- **Response prediction**: How vessels/structures will respond

MarineLabs uses spectral analysis to provide detailed wave forecasts, helping operators understand not just wave height, but the complete wave environment.

### Part 1: Understanding Wave Spectra

#### Key Concepts

**Frequency Domain Representation:**
- Wave spectrum S(f): Energy density as function of frequency
- Peak frequency fp: Frequency with maximum energy
- Spectral moments: Mn = ∫ f^n S(f) df

**Important Spectral Parameters:**
- Significant wave height: Hs = 4√(m0)
- Mean wave period: Tm02 = √(m0/m2)
- Peak period: Tp = 1/fp

### Part 2: Implementing FFT Analysis

```rust
// In src/spectral_analysis.rs
use rustfft::{FftPlanner, num_complex::Complex};
use std::f64::consts::PI;

pub struct SpectralAnalysis {
    pub frequencies: Vec<f64>,
    pub spectrum: Vec<f64>,
    pub peak_frequency: f64,
    pub peak_period: f64,
}

impl WaveTimeSeries {
    /// Perform spectral analysis using FFT
    pub fn spectral_analysis(&self) -> Result<SpectralAnalysis, WaveError> {
        if self.measurements.len() < 256 {
            return Err(WaveError::InsufficientData);
        }

        // Detrend data
        let mean = self.mean_level()?;
        let detrended: Vec<f64> = self.measurements
            .iter()
            .map(|m| m.elevation - mean)
            .collect();

        // Apply Hann window
        let n = detrended.len();
        let windowed: Vec<Complex<f64>> = detrended
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                let window = 0.5 * (1.0 - (2.0 * PI * i as f64 / (n - 1) as f64).cos());
                Complex::new(x * window, 0.0)
            })
            .collect();

        // Perform FFT
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(n);
        let mut spectrum_data = windowed.clone();
        fft.process(&mut spectrum_data);

        // Calculate power spectral density
        let df = self.sampling_rate / n as f64;
        let frequencies: Vec<f64> = (0..n/2)
            .map(|i| i as f64 * df)
            .collect();

        let spectrum: Vec<f64> = spectrum_data[..n/2]
            .iter()
            .map(|c| 2.0 * c.norm_sqr() / (self.sampling_rate * n as f64))
            .collect();

        // Find peak frequency
        let (peak_idx, _) = spectrum
            .iter()
            .enumerate()
            .skip(1) // Skip DC component
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();

        let peak_frequency = frequencies[peak_idx];
        let peak_period = 1.0 / peak_frequency;

        Ok(SpectralAnalysis {
            frequencies,
            spectrum,
            peak_frequency,
            peak_period,
        })
    }

    /// Calculate spectral moments
    pub fn spectral_moments(&self, order: i32) -> Result<f64, WaveError> {
        let analysis = self.spectral_analysis()?;
        
        let moment = analysis.frequencies
            .iter()
            .zip(analysis.spectrum.iter())
            .skip(1) // Skip DC
            .map(|(f, s)| f.powi(order) * s * (analysis.frequencies[1] - analysis.frequencies[0]))
            .sum();

        Ok(moment)
    }
}
```

### Part 3: Sea vs Swell Separation

```rust
pub struct WaveComponents {
    pub sea_spectrum: Vec<f64>,
    pub swell_spectrum: Vec<f64>,
    pub sea_hs: f64,
    pub swell_hs: f64,
    pub separation_frequency: f64,
}

impl SpectralAnalysis {
    /// Separate sea and swell components
    pub fn separate_sea_swell(&self, wind_speed: f64) -> WaveComponents {
        // Pierson-Moskowitz frequency for fully developed sea
        let fm = 0.13 * 9.81 / (2.0 * PI * wind_speed);
        
        // Separation frequency (typically 1.2 * fm)
        let separation_frequency = 1.2 * fm;
        
        // Split spectrum
        let mut sea_spectrum = vec![0.0; self.spectrum.len()];
        let mut swell_spectrum = vec![0.0; self.spectrum.len()];
        
        for (i, (f, s)) in self.frequencies.iter()
            .zip(self.spectrum.iter())
            .enumerate() {
            if *f > separation_frequency {
                sea_spectrum[i] = *s;
            } else {
                swell_spectrum[i] = *s;
            }
        }
        
        // Calculate wave heights
        let df = self.frequencies[1] - self.frequencies[0];
        let sea_m0 = sea_spectrum.iter().sum::<f64>() * df;
        let swell_m0 = swell_spectrum.iter().sum::<f64>() * df;
        
        WaveComponents {
            sea_spectrum,
            swell_spectrum,
            sea_hs: 4.0 * sea_m0.sqrt(),
            swell_hs: 4.0 * swell_m0.sqrt(),
            separation_frequency,
        }
    }
}
```

### Part 4: Directional Spectra (Advanced)

For complete wave characterization, we need directional information:

```rust
pub struct DirectionalSpectrum {
    pub frequencies: Vec<f64>,
    pub directions: Vec<f64>,    // degrees
    pub spectrum: Vec<Vec<f64>>, // S(f, θ)
}

// This requires multiple sensors or special measurement techniques
// We'll simulate for learning purposes
```

### Exercises

#### Exercise 1: JONSWAP Spectrum
Implement the JONSWAP spectrum model used for fetch-limited seas:

```rust
pub fn jonswap_spectrum(
    frequency: f64,
    hs: f64,
    peak_period: f64,
    gamma: f64,  // Peak enhancement factor (typically 3.3)
) -> f64 {
    // TODO: Implement JONSWAP formula
    // S(f) = α * g^2 / (2π)^4 / f^5 * exp(-5/4 * (fp/f)^4) * γ^exp(...)
}
```

#### Exercise 2: Wave Group Analysis
Analyze wave groupiness using spectral bandwidth:

```rust
pub struct WaveGroupParameters {
    pub groupiness_factor: f64,
    pub spectral_width: f64,
}

impl SpectralAnalysis {
    pub fn wave_group_analysis(&self) -> WaveGroupParameters {
        // TODO: Calculate spectral width parameters
        // Use spectral moments to determine groupiness
    }
}
```

#### Exercise 3: Real-time Spectral Updates
Implement a sliding window spectral analyzer:

```rust
pub struct RealTimeSpectralAnalyzer {
    window_size: usize,
    overlap: f64,
    // TODO: Add fields for streaming analysis
}

impl RealTimeSpectralAnalyzer {
    pub fn update(&mut self, new_data: &[f64]) -> Option<SpectralAnalysis> {
        // TODO: Implement sliding window FFT
        // Return updated spectrum when window is full
    }
}
```

### Questions for Reflection

1. **Vessel Safety**: How do different spectral shapes affect vessel motion differently?
