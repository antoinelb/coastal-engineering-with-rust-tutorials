use chrono::{DateTime, Datelike, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// wave observation with directional information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WaveObservation {
    pub timestamp: DateTime<Utc>,
    pub hs: f64,                 // significant wave height (m)
    pub tp: f64,                 // peak period (s)
    pub tm: f64,                 // mean period (s)
    pub dp: f64,                 // peak direction (degrees from North)
    pub dm: f64,                 // mean direction (degrees from North)
    pub spread: f64,             // directional spread (degrees)
    pub wind_speed: Option<f64>, // wind speed (m/s)
    pub wind_dir: Option<f64>,   // wind direction (degrees)
}

// wave climate statistics for a given period
#[derive(Debug, Clone)]
pub struct WaveClimateStats {
    pub mean_hs: f64,
    pub max_hs: f64,
    pub percentiles: HashMap<u8, f64>,
    pub mean_tp: f64,
    pub dominant_direction: f64,
    pub directional_spread: f64,
    pub storm_count: usize,   // number of storms where hs > storm threshold
    pub calm_percentage: f64, // percentage of time hs < calm threshold
}

// seasonal wave climate analysis
pub struct SeasonalAnalysis {
    pub winter: WaveClimateStats, // dec, jan, feb
    pub spring: WaveClimateStats, // mar, apr, may
    pub summer: WaveClimateStats, // jun, jul, aug
    pub autumn: WaveClimateStats, // sep, oct, nov
}

// directional bin for wave rose
#[derive(Debug, Clone)]
pub struct DirectionalBin {
    pub direction: f64,       // center direction (degrees)
    pub width: f64,           // bin width (degrees)
    pub hs_mean: f64,         // mean hs in this direction
    pub occurrence: f64,      // percentage of time
    pub energy_fraction: f64, // fraction of total wave energy
}

// wave climate analyzer
pub struct WaveClimateAnalyzer {
    pub observations: Vec<WaveObservation>,
    pub location: String,
    pub storm_threshold: f64,
    pub calm_threshold: f64,
}

impl WaveClimateAnalyzer {
    // create new analyzer with data
    pub fn new(location: String, storm_threshold: f64, calm_threshold: f64) -> Self {
        WaveClimateAnalyzer {
            observations: Vec::new(),
            location,
            storm_threshold,
            calm_threshold,
        }
    }

    // add observations to dataset
    pub fn add_observations(&mut self, mut obs: Vec<WaveObservation>) {
        self.observations.append(&mut obs);
        self.observations.sort_by_key(|o| o.timestamp);
    }

    // calculate basic statistics for a subset of observations
    fn calculate_stats(&self, subset: &[WaveObservation]) -> WaveClimateStats {
        if subset.is_empty() {
            panic!("Cannot calculate statistics for an empty dataset.")
        }

        // extract hs values and sort for percentiles
        let mut hs_values: Vec<f64> = subset.iter().map(|o| o.hs).collect();
        hs_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // calculate percentiles
        let mut percentiles = HashMap::new();
        for p in &[50, 90, 95, 99] {
            let idx = ((*p as f64 / 100.0) * (hs_values.len() - 1) as f64) as usize;
            percentiles.insert(*p, hs_values[idx]);
        }

        // mean values
        let mean_hs = hs_values.iter().sum::<f64>() / hs_values.len() as f64;
        let max_hs = *hs_values.last().unwrap();
        let mean_tp = subset.iter().map(|o| o.tp).sum::<f64>() / subset.len() as f64;

        // dominant direction (vector averaging)
        let (sin_sum, cos_sum) = subset.iter().fold((0.0, 0.0), |(sin_acc, cos_acc), obs| {
            let rad = obs.dp.to_radians();
            (sin_acc + rad.sin(), cos_acc + rad.cos())
        });
        let dominant_direction = sin_sum.atan2(cos_sum).to_degrees();
        let dominant_direction = if dominant_direction < 0.0 {
            dominant_direction + 360.0
        } else {
            dominant_direction
        };

        // storm and calm statistics
        let storm_count = subset
            .iter()
            .filter(|o| o.hs > self.storm_threshold)
            .count();
        let calm_count = subset.iter().filter(|o| o.hs < self.calm_threshold).count();
        let calm_percentage = (calm_count as f64 / subset.len() as f64) * 100.0;

        // average directional spread
        let directional_spread = subset.iter().map(|o| o.spread).sum::<f64>() / subset.len() as f64;

        WaveClimateStats {
            mean_hs,
            max_hs,
            percentiles,
            mean_tp,
            dominant_direction,
            directional_spread,
            storm_count,
            calm_percentage,
        }
    }

    // perform seasonal analysis
    pub fn seasonal_analysis(&self) -> SeasonalAnalysis {
        let mut winter = Vec::new();
        let mut spring = Vec::new();
        let mut summer = Vec::new();
        let mut autumn = Vec::new();

        for obs in &self.observations {
            match obs.timestamp.month() {
                12 | 1 | 2 => winter.push(obs.clone()),
                3 | 4 | 5 => spring.push(obs.clone()),
                6 | 7 | 8 => summer.push(obs.clone()),
                9 | 10 | 11 => autumn.push(obs.clone()),
                _ => unreachable!(),
            }
        }

        SeasonalAnalysis {
            winter: self.calculate_stats(&winter),
            spring: self.calculate_stats(&spring),
            summer: self.calculate_stats(&summer),
            autumn: self.calculate_stats(&autumn),
        }
    }

    // create wave rose data
    pub fn wave_rose(&self, n_bins: usize) -> Vec<DirectionalBin> {
        let bin_width = 360.0 / n_bins as f64;
        let mut bins = vec![Vec::new(); n_bins];

        // sort observations into directional bins
        for obs in &self.observations {
            let mut dir = obs.dp;
            if dir < 0.0 {
                dir += 360.0
            } else if dir >= 360.0 {
                dir -= 360.0
            };
            let bin_idx = ((dir / bin_width) as usize) % n_bins;
            bins[bin_idx].push(obs);
        }

        // calculate statistics for each bin
        let total_energy: f64 = self
            .observations
            .iter()
            .map(|o| o.hs.powi(2)) // wave energy ∝ H²
            .sum();

        bins.into_iter()
            .enumerate()
            .map(|(i, obs_vec)| {
                let direction = i as f64 * bin_width + bin_width / 2.0;

                if obs_vec.is_empty() {
                    DirectionalBin {
                        direction,
                        width: bin_width,
                        hs_mean: 0.0,
                        occurrence: 0.0,
                        energy_fraction: 0.0,
                    }
                } else {
                    let hs_mean = obs_vec.iter().map(|o| o.hs).sum::<f64>() / obs_vec.len() as f64;
                    let occurrence =
                        (obs_vec.len() as f64 / self.observations.len() as f64) * 100.0;
                    let bin_energy: f64 = obs_vec.iter().map(|o| o.hs.powi(2)).sum();
                    let energy_fraction = bin_energy / total_energy;
                    DirectionalBin {
                        direction,
                        width: bin_width,
                        hs_mean,
                        occurrence,
                        energy_fraction,
                    }
                }
            })
            .collect()
    }

    // identify swell vs wind sea components
    pub fn classify_wave_system(&self) -> (Vec<WaveObservation>, Vec<WaveObservation>) {
        let mut swell = Vec::new();
        let mut sea = Vec::new();

        for obs in &self.observations {
            // use wave age criterion: swell if tp > 1.2 * wind period
            // approximate wind period from wind speed: tw ≈ 0.33 * U (deep water)
            if let Some(wind_speed) = obs.wind_speed {
                let wind_period = 0.33 * wind_speed;
                let wave_age = obs.tp / wind_period;
                if wave_age > 1.2 {
                    swell.push(obs.clone());
                } else {
                    sea.push(obs.clone());
                }
            } else {
                // without wind data, use period threshold
                // typical BC coast: swell if tp > 10s
                if obs.tp > 10.0 {
                    swell.push(obs.clone());
                } else {
                    sea.push(obs.clone());
                }
            }
        }

        (swell, sea)
    }
}
