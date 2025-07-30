# Tutorial 4.3: Seasonal Wave Patterns
## Global Wave Environments Chapter 4.2-4.4

### Learning Objectives
By the end of this tutorial, you will:
1. Analyze seasonal variations in wave climates and their driving mechanisms
2. Implement spectral decomposition to identify wave systems (sea vs swell)
3. Quantify seasonal morphological responses to changing wave conditions
4. Model the influence of large-scale climate patterns (ENSO, PDO) on waves
5. Design coastal projects that account for seasonal variability

### Prerequisites
- Read Coastal Dynamics Ch. 4.2-4.4
- Complete Tutorials 4.1-4.2 (Wave Climate and Extreme Values)
- Understanding of atmospheric and oceanic circulation patterns
- Basic knowledge of beach morphodynamics

### Introduction: The Rhythm of the Coast

Coastal environments pulse with seasonal rhythms. Winter storms reshape beaches, while summer calms allow recovery. This seasonal dance between erosion and accretion shapes not just the physical coast but also the ecological communities and human activities that depend on it.

For the BC coast, seasonal patterns are particularly pronounced:
- **Winter (Oct-Mar)**: Intense North Pacific storms generate large waves from the W-NW
- **Summer (Apr-Sep)**: Calmer conditions with smaller waves, more from the SW
- **Transitions**: Rapid changes in wave climate during spring and fall

Understanding these patterns is crucial for:
- **Beach management**: Timing nourishment projects and assessing recovery
- **Port operations**: Scheduling maintenance during calm periods
- **Ecological protection**: Aligning activities with species life cycles
- **Indigenous practices**: Traditional knowledge of seasonal ocean conditions
- **Climate adaptation**: Detecting changes in seasonal patterns

### Part 1: Atmospheric and Oceanic Drivers

#### Global Wind Systems

The seasonal wave patterns result from large-scale atmospheric circulation:

1. **Hadley Cells**: Rising air at equator, descending at ~30°
2. **Ferrel Cells**: Mid-latitude circulation (30-60°)
3. **Polar Cells**: Cold air sinking at poles

These create the major wind belts:
- **Trade Winds**: Easterlies in tropics
- **Westerlies**: Prevailing winds at 40-60° ("Roaring Forties")
- **Polar Easterlies**: Cold winds from poles

#### Seasonal Migration of Wind Systems

The Inter-Tropical Convergence Zone (ITCZ) and storm tracks migrate seasonally:
- **Northern summer**: Systems shift northward
- **Northern winter**: Systems shift southward

This migration drives the seasonal wave patterns observed globally.

#### Climate Oscillations

Large-scale climate patterns modulate seasonal patterns:
- **ENSO** (El Niño-Southern Oscillation): 2-7 year cycle
- **PDO** (Pacific Decadal Oscillation): 20-30 year cycle
- **NAO** (North Atlantic Oscillation): Affects Atlantic coasts

#### Questions to Consider
1. How do storm tracks differ between El Niño and La Niña years?
2. Why are seasonal patterns more pronounced at higher latitudes?
3. How might climate change alter the seasonal migration of wind systems?

### Part 2: Seasonal Analysis Framework

Create a comprehensive seasonal analysis system in `src/seasonal_patterns.rs`:

```rust
// seasonal_patterns.rs
use chrono::{DateTime, Datelike, Utc, Month};
use ndarray::{Array1, Array2, Axis};
use std::collections::HashMap;
use std::f64::consts::PI;

/// Seasonal wave characteristics
#[derive(Debug, Clone)]
pub struct SeasonalCharacteristics {
    pub season: Season,
    pub mean_hs: f64,
    pub mean_tp: f64,
    pub dominant_dir: f64,
    pub directional_spread: f64,
    pub storm_frequency: f64,    // storms per month
    pub calm_frequency: f64,      // % time below threshold
    pub energy_flux: f64,         // Wave power
    pub spectral_width: f64,      // Narrow = swell, Wide = sea
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Season {
    Winter,  // Dec, Jan, Feb
    Spring,  // Mar, Apr, May
    Summer,  // Jun, Jul, Aug
    Fall,    // Sep, Oct, Nov
}

/// Monthly wave statistics
#[derive(Debug, Clone)]
pub struct MonthlyStats {
    pub month: Month,
    pub year: i32,
    pub mean_hs: f64,
    pub max_hs: f64,
    pub mean_tp: f64,
    pub dominant_dir: f64,
    pub energy_content: f64,
    pub n_storms: usize,
}

/// Seasonal pattern analyzer
pub struct SeasonalAnalyzer {
    observations: Vec<WaveObservation>,
    location: String,
    storm_threshold: f64,
    calm_threshold: f64,
}

impl SeasonalAnalyzer {
    pub fn new(location: String, storm_threshold: f64, calm_threshold: f64) -> Self {
        SeasonalAnalyzer {
            observations: Vec::new(),
            location,
            storm_threshold,
            calm_threshold,
        }
    }
    
    /// Calculate monthly statistics
    pub fn monthly_analysis(&self) -> Vec<MonthlyStats> {
        // Group by year and month
        let mut monthly_groups: HashMap<(i32, Month), Vec<&WaveObservation>> = HashMap::new();
        
        for obs in &self.observations {
            let key = (obs.timestamp.year(), Month::try_from(obs.timestamp.month() as u8).unwrap());
            monthly_groups.entry(key).or_insert_with(Vec::new).push(obs);
        }
        
        // Calculate statistics for each month
        let mut results = Vec::new();
        
        for ((year, month), obs_vec) in monthly_groups {
            if obs_vec.is_empty() { continue; }
            
            let hs_values: Vec<f64> = obs_vec.iter().map(|o| o.hs).collect();
            let mean_hs = hs_values.iter().sum::<f64>() / hs_values.len() as f64;
            let max_hs = hs_values.iter().fold(0.0, |max, &val| max.max(val));
            
            let mean_tp = obs_vec.iter().map(|o| o.tp).sum::<f64>() / obs_vec.len() as f64;
            
            // Dominant direction using vector averaging
            let (sin_sum, cos_sum) = obs_vec.iter()
                .map(|o| o.dp.to_radians())
                .fold((0.0, 0.0), |(s, c), angle| (s + angle.sin(), c + angle.cos()));
            let dominant_dir = sin_sum.atan2(cos_sum).to_degrees();
            let dominant_dir = if dominant_dir < 0.0 { dominant_dir + 360.0 } else { dominant_dir };
            
            // Wave energy (proportional to Hs²)
            let energy_content = hs_values.iter().map(|h| h * h).sum::<f64>() / hs_values.len() as f64;
            
            // Storm count
            let n_storms = hs_values.iter().filter(|&&h| h > self.storm_threshold).count();
            
            results.push(MonthlyStats {
                month,
                year,
                mean_hs,
                max_hs,
                mean_tp,
                dominant_dir,
                energy_content,
                n_storms,
            });
        }
        
        results.sort_by_key(|s| (s.year, s.month as u8));
        results
    }
    
    /// Get season for a given month
    fn get_season(month: Month) -> Season {
        match month {
            Month::December | Month::January | Month::February => Season::Winter,
            Month::March | Month::April | Month::May => Season::Spring,
            Month::June | Month::July | Month::August => Season::Summer,
            Month::September | Month::October | Month::November => Season::Fall,
        }
    }
    
    /// Calculate seasonal characteristics
    pub fn seasonal_characteristics(&self) -> HashMap<Season, SeasonalCharacteristics> {
        let mut seasonal_data: HashMap<Season, Vec<&WaveObservation>> = HashMap::new();
        
        // Group observations by season
        for obs in &self.observations {
            let month = Month::try_from(obs.timestamp.month() as u8).unwrap();
            let season = Self::get_season(month);
            seasonal_data.entry(season).or_insert_with(Vec::new).push(obs);
        }
        
        // Calculate characteristics for each season
        let mut results = HashMap::new();
        
        for (season, obs_vec) in seasonal_data {
            let mean_hs = obs_vec.iter().map(|o| o.hs).sum::<f64>() / obs_vec.len() as f64;
            let mean_tp = obs_vec.iter().map(|o| o.tp).sum::<f64>() / obs_vec.len() as f64;
            
            // Dominant direction
            let (sin_sum, cos_sum) = obs_vec.iter()
                .map(|o| o.dp.to_radians())
                .fold((0.0, 0.0), |(s, c), angle| (s + angle.sin(), c + angle.cos()));
            let dominant_dir = sin_sum.atan2(cos_sum).to_degrees();
            let dominant_dir = if dominant_dir < 0.0 { dominant_dir + 360.0 } else { dominant_dir };
            
            // Directional spread
            let directional_spread = obs_vec.iter().map(|o| o.spread).sum::<f64>() / obs_vec.len() as f64;
            
            // Storm and calm frequencies
            let n_storms = obs_vec.iter().filter(|o| o.hs > self.storm_threshold).count();
            let n_calm = obs_vec.iter().filter(|o| o.hs < self.calm_threshold).count();
            let storm_frequency = n_storms as f64 / (obs_vec.len() as f64 / (24.0 / 3.0 * 30.0)); // per month
            let calm_frequency = n_calm as f64 / obs_vec.len() as f64 * 100.0;
            
            // Wave energy flux (deep water): P = ρg²/(64π) * Hs² * Tp
            let rho = 1025.0;  // kg/m³
            let g = 9.81;      // m/s²
            let energy_flux = obs_vec.iter()
                .map(|o| rho * g * g / (64.0 * PI) * o.hs * o.hs * o.tp)
                .sum::<f64>() / obs_vec.len() as f64 / 1000.0;  // kW/m
            
            // Spectral width (using Tp variation as proxy)
            let tp_mean = mean_tp;
            let tp_std = obs_vec.iter()
                .map(|o| (o.tp - tp_mean).powi(2))
                .sum::<f64>()
                .sqrt() / obs_vec.len() as f64;
            let spectral_width = tp_std / tp_mean;
            
            results.insert(season, SeasonalCharacteristics {
                season,
                mean_hs,
                mean_tp,
                dominant_dir,
                directional_spread,
                storm_frequency,
                calm_frequency,
                energy_flux,
                spectral_width,
            });
        }
        
        results
    }
}

/// Wave observation with spectral information
#[derive(Debug, Clone)]
pub struct WaveObservation {
    pub timestamp: DateTime<Utc>,
    pub hs: f64,
    pub tp: f64,
    pub tm: f64,
    pub dp: f64,
    pub dm: f64,
    pub spread: f64,
    pub wind_speed: Option<f64>,
    pub wind_dir: Option<f64>,
}
```

### Part 3: Seasonal Morphological Response

Model how beaches respond to seasonal wave changes:

```rust
// In src/morphological_response.rs
use ndarray::{Array1, Array2};

/// Beach state classification (Wright and Short, 1984)
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BeachState {
    Dissipative,           // Ω > 6
    IntermediateLTT,       // Longshore Bar-Trough
    IntermediateRBB,       // Rhythmic Bar and Beach
    IntermediateTBR,       // Transverse Bar and Rip
    IntermediateLTRR,      // Low Tide Ridge and Runnel
    Reflective,            // Ω < 1
}

/// Dean parameter (dimensionless fall velocity)
pub fn dean_parameter(hs: f64, tp: f64, ws: f64) -> f64 {
    // Ω = Hs / (ws * T)
    // where ws is sediment fall velocity
    hs / (ws * tp)
}

/// Seasonal beach profile analyzer
pub struct BeachProfileAnalyzer {
    grain_size: f64,        // D50 in mm
    beach_slope: f64,       // Average beach face slope
    closure_depth: f64,     // Depth of closure
}

impl BeachProfileAnalyzer {
    pub fn new(grain_size: f64, beach_slope: f64, closure_depth: f64) -> Self {
        BeachProfileAnalyzer {
            grain_size,
            beach_slope,
            closure_depth,
        }
    }
    
    /// Calculate sediment fall velocity (Soulsby, 1997)
    pub fn fall_velocity(&self) -> f64 {
        let d = self.grain_size / 1000.0;  // Convert to meters
        let nu = 1.36e-6;  // Kinematic viscosity of seawater at 10°C
        let s = 2.65;      // Specific gravity of quartz sand
        let g = 9.81;
        
        // Dimensionless grain size
        let d_star = d * ((s - 1.0) * g / (nu * nu)).powf(1.0/3.0);
        
        // Fall velocity
        let ws = nu / d * (10.36_f64.sqrt() + 1.049 * d_star.powf(3.0)).sqrt() - 10.36.sqrt();
        ws
    }
    
    /// Classify beach state based on wave conditions
    pub fn classify_beach_state(&self, hs: f64, tp: f64) -> BeachState {
        let ws = self.fall_velocity();
        let omega = dean_parameter(hs, tp, ws);
        
        match omega {
            o if o < 1.0 => BeachState::Reflective,
            o if o < 2.0 => BeachState::IntermediateLTRR,
            o if o < 3.0 => BeachState::IntermediateTBR,
            o if o < 4.0 => BeachState::IntermediateRBB,
            o if o < 5.0 => BeachState::IntermediateLTT,
            _ => BeachState::Dissipative,
        }
    }
    
    /// Estimate seasonal profile change (simplified Dean equilibrium profile)
    pub fn equilibrium_profile(&self, wave_conditions: &SeasonalCharacteristics) -> Vec<(f64, f64)> {
        let ws = self.fall_velocity();
        let a = 0.067 * ws.powf(0.44);  // Dean's A parameter
        
        // Generate profile from shoreline to closure depth
        let mut profile = Vec::new();
        let x_max = (self.closure_depth / a).powf(1.5);  // Distance to closure depth
        
        for i in 0..100 {
            let x = i as f64 * x_max / 99.0;
            let h = a * x.powf(2.0/3.0);
            profile.push((x, -h));  // Negative for depth below MSL
        }
        
        profile
    }
    
    /// Calculate potential longshore transport
    pub fn longshore_transport_potential(
        &self,
        hs: f64,
        tp: f64,
        wave_angle: f64,  // Angle to shore normal in degrees
    ) -> f64 {
        // CERC formula: Q = K * Hs^2.5 * sin(2θ)
        let k = 0.39;  // CERC coefficient (SI units)
        let theta = wave_angle.to_radians();
        let g = 9.81;
        let cgb = (g * hs).sqrt();  // Approximate group velocity at breaking
        
        // Volumetric transport rate (m³/s)
        k * hs.powf(2.5) * cgb.sqrt() * (2.0 * theta).sin()
    }
    
    /// Estimate seasonal volume changes
    pub fn seasonal_volume_change(
        &self,
        summer_waves: &SeasonalCharacteristics,
        winter_waves: &SeasonalCharacteristics,
        beach_length: f64,
    ) -> VolumeChangeEstimate {
        // Summer profile (typically more reflective)
        let summer_profile = self.equilibrium_profile(summer_waves);
        let summer_area = self.profile_area(&summer_profile);
        
        // Winter profile (typically more dissipative)
        let winter_profile = self.equilibrium_profile(winter_waves);
        let winter_area = self.profile_area(&winter_profile);
        
        // Volume change per unit width
        let volume_change_per_width = summer_area - winter_area;
        let total_volume_change = volume_change_per_width * beach_length;
        
        // Estimate berm retreat/advance
        let berm_height = 2.0;  // Typical berm height above MSL
        let horizontal_change = volume_change_per_width / berm_height;
        
        VolumeChangeEstimate {
            volume_per_width: volume_change_per_width,
            total_volume: total_volume_change,
            shoreline_change: horizontal_change,
            summer_state: self.classify_beach_state(summer_waves.mean_hs, summer_waves.mean_tp),
            winter_state: self.classify_beach_state(winter_waves.mean_hs, winter_waves.mean_tp),
        }
    }
    
    /// Calculate profile area using trapezoidal rule
    fn profile_area(&self, profile: &[(f64, f64)]) -> f64 {
        let mut area = 0.0;
        for i in 1..profile.len() {
            let dx = profile[i].0 - profile[i-1].0;
            let avg_depth = (profile[i].1.abs() + profile[i-1].1.abs()) / 2.0;
            area += dx * avg_depth;
        }
        area
    }
}

#[derive(Debug)]
pub struct VolumeChangeEstimate {
    pub volume_per_width: f64,     // m³/m
    pub total_volume: f64,         // m³
    pub shoreline_change: f64,     // m (positive = accretion)
    pub summer_state: BeachState,
    pub winter_state: BeachState,
}
```

### Part 4: Climate Oscillation Analysis

Analyze the influence of ENSO and PDO on seasonal patterns:

```rust
// In src/climate_oscillations.rs
use chrono::{DateTime, Utc};
use std::collections::HashMap;

/// Climate index values
#[derive(Debug, Clone)]
pub struct ClimateIndex {
    pub timestamp: DateTime<Utc>,
    pub oni: f64,      // Oceanic Niño Index
    pub pdo: f64,      // Pacific Decadal Oscillation
    pub nao: f64,      // North Atlantic Oscillation
}

/// ENSO phase classification
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EnsoPhase {
    ElNino,    // ONI > 0.5
    Neutral,   // -0.5 ≤ ONI ≤ 0.5
    LaNina,    // ONI < -0.5
}

/// Climate pattern analyzer
pub struct ClimatePatternAnalyzer {
    climate_indices: Vec<ClimateIndex>,
    wave_observations: Vec<WaveObservation>,
}

impl ClimatePatternAnalyzer {
    /// Classify ENSO phase
    pub fn classify_enso(oni: f64) -> EnsoPhase {
        if oni > 0.5 {
            EnsoPhase::ElNino
        } else if oni < -0.5 {
            EnsoPhase::LaNina
        } else {
            EnsoPhase::Neutral
        }
    }
    
    /// Analyze wave climate by ENSO phase
    pub fn analyze_by_enso(&self) -> HashMap<EnsoPhase, WaveClimateStats> {
        let mut phase_groups: HashMap<EnsoPhase, Vec<&WaveObservation>> = HashMap::new();
        
        // Match wave observations with climate indices
        for obs in &self.wave_observations {
            // Find closest climate index (within same month)
            if let Some(index) = self.find_closest_index(&obs.timestamp) {
                let phase = Self::classify_enso(index.oni);
                phase_groups.entry(phase).or_insert_with(Vec::new).push(obs);
            }
        }
        
        // Calculate statistics for each phase
        let mut results = HashMap::new();
        
        for (phase, observations) in phase_groups {
            let stats = self.calculate_wave_stats(&observations);
            results.insert(phase, stats);
        }
        
        results
    }
    
    /// Find climate index for given date
    fn find_closest_index(&self, date: &DateTime<Utc>) -> Option<&ClimateIndex> {
        self.climate_indices.iter()
            .min_by_key(|idx| (idx.timestamp.signed_duration_since(*date)).num_days().abs())
            .filter(|idx| (idx.timestamp.signed_duration_since(*date)).num_days().abs() < 30)
    }
    
    /// Calculate composite wave statistics
    fn calculate_wave_stats(&self, observations: &[&WaveObservation]) -> WaveClimateStats {
        let mean_hs = observations.iter().map(|o| o.hs).sum::<f64>() / observations.len() as f64;
        let percentile_95 = self.percentile(observations, 95.0);
        let storm_frequency = observations.iter()
            .filter(|o| o.hs > 4.0)
            .count() as f64 / observations.len() as f64 * 100.0;
        
        // Direction statistics
        let (sin_sum, cos_sum) = observations.iter()
            .map(|o| o.dp.to_radians())
            .fold((0.0, 0.0), |(s, c), angle| (s + angle.sin(), c + angle.cos()));
        let mean_direction = sin_sum.atan2(cos_sum).to_degrees();
        
        WaveClimateStats {
            mean_hs,
            percentile_95,
            storm_frequency,
            mean_direction: if mean_direction < 0.0 { mean_direction + 360.0 } else { mean_direction },
            n_observations: observations.len(),
        }
    }
    
    fn percentile(&self, observations: &[&WaveObservation], p: f64) -> f64 {
        let mut hs_values: Vec<f64> = observations.iter().map(|o| o.hs).collect();
        hs_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let idx = ((p / 100.0) * (hs_values.len() - 1) as f64) as usize;
        hs_values[idx]
    }
}

#[derive(Debug)]
pub struct WaveClimateStats {
    pub mean_hs: f64,
    pub percentile_95: f64,
    pub storm_frequency: f64,
    pub mean_direction: f64,
    pub n_observations: usize,
}

/// Teleconnection pattern analyzer
pub struct TeleconnectionAnalyzer;

impl TeleconnectionAnalyzer {
    /// Calculate correlation between climate index and wave parameters
    pub fn calculate_correlation(
        index_values: &[f64],
        wave_parameter: &[f64],
        lag_months: i32,
    ) -> f64 {
        // Apply lag
        let (index_lagged, wave_aligned) = if lag_months >= 0 {
            let lag = lag_months as usize;
            (&index_values[..index_values.len()-lag], &wave_parameter[lag..])
        } else {
            let lag = (-lag_months) as usize;
            (&index_values[lag..], &wave_parameter[..wave_parameter.len()-lag])
        };
        
        // Calculate Pearson correlation
        let n = index_lagged.len().min(wave_aligned.len()) as f64;
        let mean_x = index_lagged.iter().sum::<f64>() / n;
        let mean_y = wave_aligned.iter().sum::<f64>() / n;
        
        let cov = index_lagged.iter().zip(wave_aligned)
            .map(|(x, y)| (x - mean_x) * (y - mean_y))
            .sum::<f64>() / n;
        
        let std_x = (index_lagged.iter()
            .map(|x| (x - mean_x).powi(2))
            .sum::<f64>() / n).sqrt();
        
        let std_y = (wave_aligned.iter()
            .map(|y| (y - mean_y).powi(2))
            .sum::<f64>() / n).sqrt();
        
        cov / (std_x * std_y)
    }
}
```

### Part 5: Visualization Functions

Create visualizations for seasonal patterns:

```rust
// In src/seasonal_visualization.rs
use plotters::prelude::*;
use crate::seasonal_patterns::*;
use crate::morphological_response::*;

/// Plot monthly wave climate
pub fn plot_monthly_climatology(
    monthly_stats: &[MonthlyStats],
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Monthly Wave Climatology", ("sans-serif", 30))
        .margin(15)
        .set_label_area_size(LabelAreaPosition::Left, 50)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(1u32..13u32, 0.0..6.0)?;
    
    chart.configure_mesh()
        .x_desc("Month")
        .y_desc("Significant Wave Height (m)")
        .x_labels(12)
        .x_label_formatter(&|x| {
            match x {
                1 => "Jan", 2 => "Feb", 3 => "Mar", 4 => "Apr",
                5 => "May", 6 => "Jun", 7 => "Jul", 8 => "Aug",
                9 => "Sep", 10 => "Oct", 11 => "Nov", 12 => "Dec",
                _ => "",
            }.to_string()
        })
        .draw()?;
    
    // Calculate monthly averages across all years
    let mut monthly_means: HashMap<u8, Vec<f64>> = HashMap::new();
    let mut monthly_maxima: HashMap<u8, Vec<f64>> = HashMap::new();
    
    for stat in monthly_stats {
        let month_num = stat.month as u8;
        monthly_means.entry(month_num).or_insert_with(Vec::new).push(stat.mean_hs);
        monthly_maxima.entry(month_num).or_insert_with(Vec::new).push(stat.max_hs);
    }
    
    // Plot mean values with error bars
    let mean_series: Vec<_> = (1..=12).map(|month| {
        let values = monthly_means.get(&month).unwrap();
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let std = (values.iter()
            .map(|v| (v - mean).powi(2))
            .sum::<f64>() / values.len() as f64).sqrt();
        (month as u32, mean, std)
    }).collect();
    
    // Plot mean line
    chart.draw_series(LineSeries::new(
        mean_series.iter().map(|(m, mean, _)| (*m, *mean)),
        &BLUE,
    ))?
    .label("Mean Hs")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &BLUE));
    
    // Add error bars (±1 std)
    chart.draw_series(
        mean_series.iter().map(|(month, mean, std)| {
            ErrorBar::new_vertical(
                (*month, mean - std),
                (*month, *mean),
                (*month, mean + std),
            )
        })
    )?;
    
    // Plot maximum envelope
    let max_series: Vec<_> = (1..=12).map(|month| {
        let values = monthly_maxima.get(&month).unwrap();
        let max = values.iter().fold(0.0, |a, &b| a.max(b));
        (month as u32, max)
    }).collect();
    
    chart.draw_series(LineSeries::new(
        max_series,
        ShapeStyle::from(&RED).stroke_width(1),
    ))?
    .label("Monthly maxima")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &RED));
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
    
    root.present()?;
    Ok(())
}

/// Plot seasonal beach profiles
pub fn plot_seasonal_profiles(
    analyzer: &BeachProfileAnalyzer,
    seasonal_chars: &HashMap<Season, SeasonalCharacteristics>,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (1000, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let summer = seasonal_chars.get(&Season::Summer).unwrap();
    let winter = seasonal_chars.get(&Season::Winter).unwrap();
    
    let summer_profile = analyzer.equilibrium_profile(summer);
    let winter_profile = analyzer.equilibrium_profile(winter);
    
    // Find plot bounds
    let x_max = summer_profile.last().unwrap().0.max(winter_profile.last().unwrap().0);
    let y_min = summer_profile.iter().map(|(_, y)| *y).fold(0.0, f64::min);
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Seasonal Beach Profile Changes", ("sans-serif", 25))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..x_max * 1.1, y_min * 1.1..2.0)?;
    
    chart.configure_mesh()
        .x_desc("Distance Offshore (m)")
        .y_desc("Elevation (m MSL)")
        .draw()?;
    
    // Add MSL line
    chart.draw_series(LineSeries::new(
        vec![(0.0, 0.0), (x_max * 1.1, 0.0)],
        ShapeStyle::from(&BLACK).stroke_width(1),
    ))?;
    
    // Plot profiles
    chart.draw_series(LineSeries::new(
        summer_profile,
        &GREEN,
    ))?
    .label(format!("Summer (Hs={:.1}m)", summer.mean_hs))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &GREEN));
    
    chart.draw_series(LineSeries::new(
        winter_profile,
        &BLUE,
    ))?
    .label(format!("Winter (Hs={:.1}m)", winter.mean_hs))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &BLUE));
    
    // Shade difference
    let mut difference_polygon = vec![];
    for (s_point, w_point) in summer_profile.iter().zip(&winter_profile) {
        difference_polygon.push(*s_point);
    }
    for w_point in winter_profile.iter().rev() {
        difference_polygon.push(*w_point);
    }
    
    chart.draw_series(std::iter::once(Polygon::new(
        difference_polygon,
        RED.mix(0.2).filled(),
    )))?;
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
    
    root.present()?;
    Ok(())
}

/// Plot ENSO influence on wave climate
pub fn plot_enso_composite(
    enso_stats: &HashMap<EnsoPhase, WaveClimateStats>,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let phases = vec![EnsoPhase::ElNino, EnsoPhase::Neutral, EnsoPhase::LaNina];
    let labels = vec!["El Niño", "Neutral", "La Niña"];
    
    let mut chart = ChartBuilder::on(&root)
        .caption("ENSO Influence on Wave Climate", ("sans-serif", 25))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(
            labels.as_slice(),
            0.0..enso_stats.values().map(|s| s.mean_hs).fold(0.0, f64::max) * 1.2,
        )?;
    
    chart.configure_mesh()
        .y_desc("Mean Significant Wave Height (m)")
        .draw()?;
    
    // Plot bars
    chart.draw_series(
        phases.iter().zip(&labels).map(|(phase, label)| {
            let stats = enso_stats.get(phase).unwrap();
            Rectangle::new([(*label, 0.0), (*label, stats.mean_hs)], BLUE.filled())
        })
    )?;
    
    // Add 95th percentile markers
    chart.draw_series(
        phases.iter().zip(&labels).map(|(phase, label)| {
            let stats = enso_stats.get(phase).unwrap();
            Circle::new((*label, stats.percentile_95), 5, RED.filled())
        })
    )?
    .label("95th percentile")
    .legend(|(x, y)| Circle::new((x + 5, y + 5), 5, RED.filled()));
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
    
    root.present()?;
    Ok(())
}
```

### Part 6: Main Application

Integrate all components:

```rust
// In src/main.rs
mod seasonal_patterns;
mod morphological_response;
mod climate_oscillations;
mod seasonal_visualization;

use seasonal_patterns::*;
use morphological_response::*;
use climate_oscillations::*;
use std::error::Error;

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {
    println!("Coastal Dynamics Tutorial 4.3: Seasonal Wave Patterns");
    println!("===================================================\n");
    
    // Load or generate multi-year wave data
    let wave_data = generate_multi_year_data(10)?;  // 10 years
    
    // Initialize seasonal analyzer
    let mut analyzer = SeasonalAnalyzer::new(
        "Tofino, BC".to_string(),
        4.0,  // Storm threshold
        0.5,  // Calm threshold
    );
    analyzer.observations = wave_data.clone();
    
    // Monthly analysis
    let monthly_stats = analyzer.monthly_analysis();
    println!("Analyzed {} months of data\n", monthly_stats.len());
    
    // Seasonal characteristics
    let seasonal_chars = analyzer.seasonal_characteristics();
    
    println!("=== Seasonal Wave Characteristics ===");
    for season in &[Season::Winter, Season::Spring, Season::Summer, Season::Fall] {
        if let Some(chars) = seasonal_chars.get(season) {
            println!("{:?}:", season);
            println!("  Mean Hs: {:.2} m", chars.mean_hs);
            println!("  Mean Tp: {:.1} s", chars.mean_tp);
            println!("  Dominant direction: {:.0}°", chars.dominant_dir);
            println!("  Storm frequency: {:.1} per month", chars.storm_frequency);
            println!("  Calm conditions: {:.1}%", chars.calm_frequency);
            println!("  Wave power: {:.1} kW/m", chars.energy_flux);
            println!("  Spectral width: {:.2} ({})",
                chars.spectral_width,
                if chars.spectral_width < 0.15 { "swell-dominated" } else { "mixed sea" }
            );
            println!();
        }
    }
    
    // Beach morphological response
    println!("=== Beach Morphological Response ===");
    let beach_analyzer = BeachProfileAnalyzer::new(
        0.25,   // D50 = 0.25 mm (fine sand)
        0.1,    // Beach slope = 1:10
        8.0,    // Closure depth = 8 m
    );
    
    let summer = seasonal_chars.get(&Season::Summer).unwrap();
    let winter = seasonal_chars.get(&Season::Winter).unwrap();
    
    let volume_change = beach_analyzer.seasonal_volume_change(
        summer,
        winter,
        1000.0,  // 1 km beach length
    );
    
    println!("Beach State Changes:");
    println!("  Summer: {:?}", volume_change.summer_state);
    println!("  Winter: {:?}", volume_change.winter_state);
    println!("\nVolume Changes:");
    println!("  Per unit width: {:.1} m³/m", volume_change.volume_per_width);
    println!("  Total (1 km beach): {:.0} m³", volume_change.total_volume);
    println!("  Shoreline change: {:.1} m", volume_change.shoreline_change);
    
    // Climate oscillation analysis
    println!("\n=== Climate Oscillation Analysis ===");
    let climate_data = load_climate_indices()?;
    
    let mut climate_analyzer = ClimatePatternAnalyzer {
        climate_indices: climate_data,
        wave_observations: wave_data,
    };
    
    let enso_stats = climate_analyzer.analyze_by_enso();
    
    for phase in &[EnsoPhase::ElNino, EnsoPhase::Neutral, EnsoPhase::LaNina] {
        if let Some(stats) = enso_stats.get(phase) {
            println!("{:?}:", phase);
            println!("  Mean Hs: {:.2} m", stats.mean_hs);
            println!("  95th percentile: {:.2} m", stats.percentile_95);
            println!("  Storm frequency: {:.1}%", stats.storm_frequency);
            println!("  Mean direction: {:.0}°", stats.mean_direction);
            println!("  Sample size: {} observations", stats.n_observations);
        }
    }
    
    // Generate visualizations
    seasonal_visualization::plot_monthly_climatology(&monthly_stats, "monthly_climate.png")?;
    seasonal_visualization::plot_seasonal_profiles(&beach_analyzer, &seasonal_chars, "seasonal_profiles.png")?;
    seasonal_visualization::plot_enso_composite(&enso_stats, "enso_composite.png")?;
    
    // Environmental and social implications
    println!("\n=== Environmental and Social Context ===");
    println!("Seasonal Ecological Cycles:");
    println!("- Winter storms create habitat for specialized species");
    println!("- Summer calms allow kelp forest recovery and growth");
    println!("- Intertidal zonation shifts with seasonal wave exposure");
    println!("- Shorebird nesting aligns with calm summer conditions");
    
    println!("\nIndigenous Knowledge:");
    println!("- Traditional calendars reflect seasonal ocean conditions");
    println!("- Harvesting practices adapt to wave exposure");
    println!("- Oral histories document extreme seasonal events");
    
    println!("\nClimate Change Implications:");
    println!("- Intensifying winter storms may exceed beach recovery capacity");
    println!("- Shifting storm tracks could alter seasonal patterns");
    println!("- Changes in ENSO frequency affect inter-annual variability");
    println!("- Coastal communities must adapt to new seasonal rhythms");
    
    Ok(())
}

/// Generate synthetic multi-year wave data with realistic patterns
fn generate_multi_year_data(years: usize) -> Result<Vec<WaveObservation>, Box<dyn Error>> {
    use rand::{thread_rng, Rng};
    use rand_distr::{Normal, Distribution};
    use chrono::{TimeZone, Duration};
    
    let mut rng = thread_rng();
    let mut observations = Vec::new();
    
    let start = chrono::Utc.ymd(2015, 1, 1).and_hms(0, 0, 0);
    
    // Generate ENSO-like pattern
    let enso_period = 5.0;  // years
    
    for day in 0..(365 * years) {
        let timestamp = start + Duration::days(day as i64);
        let day_of_year = timestamp.ordinal();
        
        // ENSO influence
        let enso_phase = (2.0 * PI * day as f64 / (365.0 * enso_period)).sin();
        
        // Seasonal pattern
        let seasonal = 1.0 + 0.8 * (2.0 * PI * (day_of_year as f64 - 80.0) / 365.0).cos();
        
        // Base conditions with ENSO modulation
        let base_hs = 1.5 * seasonal * (1.0 + 0.2 * enso_phase);
        
        // Generate 8 observations per day
        for hour in (0..24).step_by(3) {
            let timestamp = timestamp + Duration::hours(hour);
            
            let hs = Normal::new(base_hs, 0.4 * seasonal).unwrap().sample(&mut rng).max(0.1);
            let tp = 8.0 + 2.0 * hs.sqrt() + rng.gen::<f64>() * 2.0;
            let tm = tp * (0.7 + rng.gen::<f64>() * 0.2);
            
            // Seasonal direction shift
            let base_dir = if day_of_year < 90 || day_of_year > 270 {
                285.0  // Winter: WNW
            } else {
                255.0  // Summer: WSW
            };
            
            let dp = base_dir + Normal::new(0.0, 15.0).unwrap().sample(&mut rng);
            let dm = dp + Normal::new(0.0, 5.0).unwrap().sample(&mut rng);
            
            observations.push(WaveObservation {
                timestamp,
                hs,
                tp,
                tm,
                dp: dp % 360.0,
                dm: dm % 360.0,
                spread: 25.0 + rng.gen::<f64>() * 10.0,
                wind_speed: Some(5.0 + 3.0 * hs),
                wind_dir: Some(dp - 10.0),
            });
        }
    }
    
    Ok(observations)
}

/// Load climate indices (placeholder)
fn load_climate_indices() -> Result<Vec<ClimateIndex>, Box<dyn Error>> {
    // In practice, load from NOAA or similar
    // Generate synthetic data for demonstration
    use chrono::{TimeZone, Duration};
    let mut indices = Vec::new();
    
    let start = chrono::Utc.ymd(2015, 1, 1).and_hms(0, 0, 0);
    
    for month in 0..120 {  // 10 years
        let timestamp = start + Duration::days(month * 30);
        
        // Synthetic ENSO pattern
        let oni = (2.0 * PI * month as f64 / 60.0).sin() * 1.5;  // 5-year cycle
        let pdo = (2.0 * PI * month as f64 / 240.0).sin() * 1.0;  // 20-year cycle
        let nao = rand::random::<f64>() * 2.0 - 1.0;  // Random
        
        indices.push(ClimateIndex {
            timestamp,
            oni,
            pdo,
            nao,
        });
    }
    
    Ok(indices)
}
```

### Exercises

#### Exercise 1: Seasonal Sediment Budget
Implement a full seasonal sediment budget calculation:

```rust
pub struct SeasonalSedimentBudget {
    pub summer_input: f64,     // m³
    pub summer_output: f64,
    pub winter_input: f64,
    pub winter_output: f64,
    pub annual_net: f64,
}

impl BeachProfileAnalyzer {
    pub fn calculate_sediment_budget(
        &self,
        seasonal_waves: &HashMap<Season, SeasonalCharacteristics>,
        beach_geometry: &BeachGeometry,
    ) -> SeasonalSedimentBudget {
        // TODO: Calculate sources and sinks
        // - Longshore transport gradients
        // - Cross-shore exchanges
        // - Aeolian transport
        // - Cliff/bluff erosion
    }
}
```

#### Exercise 2: Spectral Partitioning
Separate sea and swell components:

```rust
pub fn partition_spectrum(
    spectrum: &WaveSpectrum,
) -> (WaveSpectrum, WaveSpectrum) {
    // TODO: Implement spectral partitioning
    // - Identify spectral peaks
    // - Separate wind sea from swell
    // - Calculate individual Hs, Tp, Dp
}
```

#### Exercise 3: Non-stationary Seasonal Analysis
Detect trends in seasonal patterns:

```rust
pub struct SeasonalTrend {
    pub hs_trend: f64,         // m/decade
    pub tp_trend: f64,         // s/decade
    pub storm_freq_trend: f64, // storms/month/decade
    pub significance: f64,     // p-value
}

pub fn analyze_seasonal_trends(
    long_term_data: &[WaveObservation],
) -> HashMap<Season, SeasonalTrend> {
    // TODO: Implement trend detection
    // - Linear regression by season
    // - Mann-Kendall test
    // - Account for climate oscillations
}
```

#### Exercise 4: Ecological Habitat Modeling
Link seasonal wave patterns to species distribution:

```rust
pub struct HabitatSuitability {
    pub species: String,
    pub seasonal_suitability: HashMap<Season, f64>,
    pub limiting_factors: Vec<String>,
}

pub fn model_habitat_suitability(
    species: &Species,
    seasonal_waves: &HashMap<Season, SeasonalCharacteristics>,
) -> HabitatSuitability {
    // TODO: Calculate habitat suitability
    // - Wave exposure tolerance
    // - Seasonal life cycle needs
    // - Food availability
}
```

### Questions for Reflection

1. **Climate Change**: How might the intensification of storms affect the balance between winter erosion and summer recovery? What are the thresholds for irreversible change?

2. **Management Timing**: When is the optimal time for beach nourishment considering seasonal patterns, ecological cycles, and human use?

3. **Traditional Knowledge**: How do Indigenous seasonal calendars align with scientific observations of ocean conditions? What can we learn from this convergence?

4. **Teleconnections**: How do remote climate patterns (ENSO, PDO) manifest in local seasonal wave conditions? What are the lag times?

5. **Morphological Memory**: How long does it take beaches to adjust to changing seasonal patterns? What determines recovery time?

### Additional Resources

- Wright and Short (1984) - Beach morphodynamics
- Ruggiero et al. (2010) - Seasonal shoreline variability
- Barnard et al. (2015) - Coastal vulnerability to climate variability
- NOAA Climate Prediction Center - ENSO resources
- Traditional Ecological Knowledge studies of Pacific Northwest

### Summary

This tutorial has explored the seasonal rhythms that shape our coasts. We've seen how global atmospheric patterns drive local wave climates, how beaches respond morphologically to seasonal changes, and how climate oscillations modulate these patterns over longer timescales.

The seasonal dance between erosion and accretion is fundamental to coastal resilience. Understanding these patterns - both through scientific measurement and traditional knowledge - is essential for sustainable coastal management in a changing climate.

Remember: The coast is not static but constantly adjusting to changing conditions. Our role as coastal engineers is to work with these natural rhythms, not against them.
