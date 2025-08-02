use crate::wave_climate::WaveObservation;
use ndarray::Array1;

#[derive(Debug)]
pub struct GevParameters {
    pub location: f64, // μ
    pub scale: f64,    // σ
    pub shape: f64,    // ξ
}

#[derive(Debug)]
pub struct PotResults {
    pub threshold: f64,
    pub rate_per_year: f64,
    pub scale: f64,
    pub n_exceedances: usize,
}

// extreme value analysis for wave heights
pub struct ExtremeValueAnalysis {
    annual_maxima: Vec<f64>,
    threshold: f64, // for POT analysis
}

impl ExtremeValueAnalysis {
    // fit generalized extreme value (GEV) distribution
    pub fn fit_gev(&self) -> GevParameters {
        // simplified GEV fitting using method of moments
        let data = Array1::from_vec(self.annual_maxima.clone());
        let mean = data.mean().unwrap();
        let std = data.std(0.0);

        // initial estimates (simplified)
        let shape = -0.1; // typical for wave heights
        let scale = std * 0.78; // approximation
        let location = mean - scale * 0.577; // Euler's constant

        GevParameters {
            location,
            scale,
            shape,
        }
    }

    // estimate return values
    pub fn return_value(&self, return_period: f64, params: &GevParameters) -> f64 {
        let p = 1.0 - 1.0 / return_period;

        if params.shape.abs() < 1e-10 {
            // Gumbel distribution (shape ≈ 0)
            params.location - params.scale * (-(-p.ln()).ln())
        } else {
            // general GEV
            params.location + params.scale / params.shape * (1.0 - (-p.ln()).powf(-params.shape))
        }
    }

    // peak-over-threshold analysis
    pub fn pot_analysis(&self, observations: &[WaveObservation]) -> PotResults {
        let exceedances: Vec<f64> = observations
            .iter()
            .filter_map(|o| {
                if o.hs > self.threshold {
                    Some(o.hs - self.threshold)
                } else {
                    None
                }
            })
            .collect();

        let n_years = observations.len() as f64 / (365.25 * 24.0 / 3.0); // assuming 3-hour data
        let rate = exceedances.len() as f64 / n_years;

        // fit exponential distribution to exceedances (simplified GPD with shape=0)
        let mean_excess = exceedances.iter().sum::<f64>() / exceedances.len() as f64;

        PotResults {
            threshold: self.threshold,
            rate_per_year: rate,
            scale: mean_excess,
            n_exceedances: exceedances.len(),
        }
    }
}
