# Tutorial 3.4: Tide Prediction Engine
## Ocean Waves Chapter 3.7-3.9

### Learning Objectives
By the end of this tutorial, you will:
1. Understand tidal generation and constituents
2. Implement harmonic analysis for tide prediction
3. Build a complete tide prediction engine
4. Apply to Canadian tidal stations
5. Integrate tides with other coastal processes

### Prerequisites
- Completed Tutorials 3.1-3.3
- Read Coastal Dynamics Ch. 3.7-3.9
- Basic understanding of harmonic analysis

### Introduction: The Rhythm of the Coast

Tides are the heartbeat of coastal systems, driving:
- **Navigation timing**: When vessels can enter/leave ports
- **Intertidal ecology**: Exposure cycles for marine life
- **Sediment transport**: Tidal currents move material
- **Renewable energy**: Predictable tidal power
- **Cultural practices**: Traditional harvesting windows

MarineLabs integrates tidal predictions with wave forecasts to provide complete water level information, crucial for safe port operations and coastal planning.

### Part 1: Understanding Tidal Forces

#### Key Concepts

**Tide-Generating Forces:**
- Gravitational attraction (Moon and Sun)
- Centrifugal force from Earth-Moon rotation
- Results in differential force across Earth

**Tidal Characteristics:**
- **Semi-diurnal**: Two high/low tides per day (Atlantic Canada)
- **Mixed**: Unequal highs/lows (Pacific Canada)
- **Diurnal**: One high/low per day (rare in Canada)

**Important Periods:**
- M2 (principal lunar): 12.42 hours
- S2 (principal solar): 12.00 hours
- K1 (diurnal): 23.93 hours
- O1 (diurnal): 25.82 hours

### Part 2: Core Tide Prediction Engine

```rust
// In src/tide_prediction.rs
use std::collections::HashMap;
use chrono::{DateTime, Utc, Duration};
use std::f64::consts::PI;

/// Tidal constituent parameters
#[derive(Debug, Clone)]
pub struct TidalConstituent {
    pub name: String,
    pub frequency: f64,      // cycles per hour
    pub amplitude: f64,      // meters
    pub phase: f64,         // degrees
    pub speed: f64,         // degrees per hour
}

/// Collection of constituents for a location
#[derive(Debug, Clone)]
pub struct TidalStation {
    pub name: String,
    pub latitude: f64,
    pub longitude: f64,
    pub timezone: i32,
    pub datum_offset: f64,   // Mean sea level to chart datum
    pub constituents: Vec<TidalConstituent>,
}

/// Tide prediction engine
pub struct TidePredictionEngine {
    stations: HashMap<String, TidalStation>,
}

impl TidePredictionEngine {
    pub fn new() -> Self {
        TidePredictionEngine {
            stations: HashMap::new(),
        }
    }
    
    /// Add a tidal station
    pub fn add_station(&mut self, station: TidalStation) {
        self.stations.insert(station.name.clone(), station);
    }
    
    /// Predict tide height at specific time
    pub fn predict_height(
        &self,
        station_name: &str,
        time: DateTime<Utc>,
    ) -> Result<f64, String> {
        let station = self.stations.get(station_name)
            .ok_or_else(|| format!("Station {} not found", station_name))?;
        
        // Reference time for phase calculation (usually start of year)
        let reference_time = DateTime::parse_from_rfc3339("2024-01-01T00:00:00Z")
            .unwrap()
            .with_timezone(&Utc);
        
        let hours_since_reference = (time - reference_time).num_seconds() as f64 / 3600.0;
        
        // Sum all constituents
        let mut height = station.datum_offset;
        
        for constituent in &station.constituents {
            let phase_rad = constituent.phase * PI / 180.0;
            let omega = constituent.speed * PI / 180.0; // Convert to radians per hour
            
            height += constituent.amplitude * 
                (omega * hours_since_reference + phase_rad).cos();
        }
        
        Ok(height)
    }
    
    /// Find high and low tides in a time range
    pub fn find_extremes(
        &self,
        station_name: &str,
        start_time: DateTime<Utc>,
        end_time: DateTime<Utc>,
    ) -> Result<Vec<TideExtreme>, String> {
        let mut extremes = Vec::new();
        let mut current_time = start_time;
        let time_step = Duration::minutes(6); // 6-minute intervals
        
        let mut previous_height = self.predict_height(station_name, current_time)?;
        current_time = current_time + time_step;
        let mut current_height = self.predict_height(station_name, current_time)?;
        let mut increasing = current_height > previous_height;
        
        while current_time < end_time {
            let next_time = current_time + time_step;
            let next_height = self.predict_height(station_name, next_time)?;
            
            let now_increasing = next_height > current_height;
            
            // Detect extremum
            if increasing != now_increasing {
                // Refine extremum time using quadratic interpolation
                let refined = self.refine_extremum(
                    station_name,
                    current_time - time_step,
                    current_time,
                    current_time + time_step,
                )?;
                
                extremes.push(TideExtreme {
                    time: refined.0,
                    height: refined.1,
                    extreme_type: if increasing { 
                        ExtremeType::High 
                    } else { 
                        ExtremeType::Low 
                    },
                });
                
                increasing = now_increasing;
            }
            
            previous_height = current_height;
            current_height = next_height;
            current_time = next_time;
        }
        
        Ok(extremes)
    }
    
    /// Refine extremum using quadratic fit
    fn refine_extremum(
        &self,
        station_name: &str,
        t1: DateTime<Utc>,
        t2: DateTime<Utc>,
        t3: DateTime<Utc>,
    ) -> Result<(DateTime<Utc>, f64), String> {
        let h1 = self.predict_height(station_name, t1)?;
        let h2 = self.predict_height(station_name, t2)?;
        let h3 = self.predict_height(station_name, t3)?;
        
        // Quadratic interpolation
        let dt = (t2 - t1).num_seconds() as f64;
        let a = (h3 - 2.0 * h2 + h1) / (2.0 * dt * dt);
        let b = (h3 - h1) / (2.0 * dt);
        
        let t_offset = -b / (2.0 * a);
        let extremum_time = t2 + Duration::seconds(t_offset as i64);
        let extremum_height = self.predict_height(station_name, extremum_time)?;
        
        Ok((extremum_time, extremum_height))
    }
}

#[derive(Debug, Clone)]
pub struct TideExtreme {
    pub time: DateTime<Utc>,
    pub height: f64,
    pub extreme_type: ExtremeType,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ExtremeType {
    High,
    Low,
}
```

### Part 3: Harmonic Analysis Implementation

```rust
/// Harmonic analysis to extract constituents from water level data
pub struct HarmonicAnalysis {
    /// Known tidal frequencies (cycles per hour)
    constituent_frequencies: HashMap<String, f64>,
}

impl HarmonicAnalysis {
    pub fn new() -> Self {
        let mut frequencies = HashMap::new();
        
        // Major constituents
        frequencies.insert("M2".to_string(), 1.0 / 12.4206012);  // Principal lunar
        frequencies.insert("S2".to_string(), 1.0 / 12.0);        // Principal solar
        frequencies.insert("N2".to_string(), 1.0 / 12.65834751); // Lunar elliptic
        frequencies.insert("K2".to_string(), 1.0 / 11.96723606); // Luni-solar
        frequencies.insert("K1".to_string(), 1.0 / 23.93447213); // Luni-solar diurnal
        frequencies.insert("O1".to_string(), 1.0 / 25.81933871); // Principal lunar diurnal
        frequencies.insert("P1".to_string(), 1.0 / 24.06588766); // Principal solar diurnal
        frequencies.insert("Q1".to_string(), 1.0 / 26.86835775); // Lunar elliptic diurnal
        
        HarmonicAnalysis {
            constituent_frequencies: frequencies,
        }
    }
    
    /// Analyze water level time series to extract constituents
    pub fn analyze(
        &self,
        water_levels: &[(DateTime<Utc>, f64)],
        constituents_to_fit: &[String],
    ) -> Result<Vec<TidalConstituent>, String> {
        if water_levels.len() < 720 {  // At least 30 days of hourly data
            return Err("Insufficient data for harmonic analysis".to_string());
        }
        
        let n = water_levels.len();
        let start_time = water_levels[0].0;
        
        // Calculate mean sea level
        let msl: f64 = water_levels.iter().map(|(_, h)| h).sum::<f64>() / n as f64;
        
        let mut constituents = Vec::new();
        
        for name in constituents_to_fit {
            let frequency = self.constituent_frequencies.get(name)
                .ok_or_else(|| format!("Unknown constituent: {}", name))?;
            
            let omega = 2.0 * PI * frequency; // radians per hour
            
            // Least squares fit: h(t) = A*cos(ωt) + B*sin(ωt)
            let mut sum_cos = 0.0;
            let mut sum_sin = 0.0;
            let mut sum_h_cos = 0.0;
            let mut sum_h_sin = 0.0;
            
            for (time, height) in water_levels {
                let t = (*time - start_time).num_seconds() as f64 / 3600.0;
                let h = height - msl;
                
                let cos_wt = (omega * t).cos();
                let sin_wt = (omega * t).sin();
                
                sum_cos += cos_wt * cos_wt;
                sum_sin += sin_wt * sin_wt;
                sum_h_cos += h * cos_wt;
                sum_h_sin += h * sin_wt;
            }
            
            let a = sum_h_cos / sum_cos;
            let b = sum_h_sin / sum_sin;
            
            // Convert to amplitude and phase
            let amplitude = (a * a + b * b).sqrt();
            let phase = b.atan2(a) * 180.0 / PI;
            
            constituents.push(TidalConstituent {
                name: name.clone(),
                frequency: *frequency,
                amplitude,
                phase,
                speed: frequency * 360.0, // degrees per hour
            });
        }
        
        Ok(constituents)
    }
}
```

### Part 4: Canadian Tidal Stations

```rust
/// Pre-defined Canadian tidal stations
pub fn load_canadian_stations() -> HashMap<String, TidalStation> {
    let mut stations = HashMap::new();
    
    // Vancouver (mixed, mainly semi-diurnal)
    stations.insert("Vancouver".to_string(), TidalStation {
        name: "Vancouver".to_string(),
        latitude: 49.2827,
        longitude: -123.1207,
        timezone: -8,
        datum_offset: 3.0, // Chart datum to MSL
        constituents: vec![
            TidalConstituent {
                name: "M2".to_string(),
                frequency: 1.0 / 12.4206012,
                amplitude: 1.20,
                phase: 235.0,
                speed: 28.984104,
            },
            TidalConstituent {
                name: "K1".to_string(),
                frequency: 1.0 / 23.93447213,
                amplitude: 0.80,
                phase: 285.0,
                speed: 15.041069,
            },
            TidalConstituent {
                name: "O1".to_string(),
                frequency: 1.0 / 25.81933871,
                amplitude: 0.50,
                phase: 265.0,
                speed: 13.943035,
            },
            TidalConstituent {
                name: "S2".to_string(),
                frequency: 1.0 / 12.0,
                amplitude: 0.35,
                phase: 270.0,
                speed: 30.0,
            },
        ],
    });
    
    // Quebec City (semi-diurnal, large range)
    stations.insert("Quebec City".to_string(), TidalStation {
        name: "Quebec City".to_string(),
        latitude: 46.8139,
        longitude: -71.2082,
        timezone: -5,
        datum_offset: 2.5,
        constituents: vec![
            TidalConstituent {
                name: "M2".to_string(),
                frequency: 1.0 / 12.4206012,
                amplitude: 2.10,
                phase: 95.0,
                speed: 28.984104,
            },
            TidalConstituent {
                name: "S2".to_string(),
                frequency: 1.0 / 12.0,
                amplitude: 0.65,
                phase: 130.0,
                speed: 30.0,
            },
            TidalConstituent {
                name: "N2".to_string(),
                frequency: 1.0 / 12.65834751,
                amplitude: 0.45,
                phase: 80.0,
                speed: 28.439730,
            },
            TidalConstituent {
                name: "K1".to_string(),
                frequency: 1.0 / 23.93447213,
                amplitude: 0.15,
                phase: 290.0,
                speed: 15.041069,
            },
        ],
    });
    
    // Add more stations as needed...
    
    stations
}
```

### Part 5: Tidal Current Prediction

```rust
/// Tidal current prediction
#[derive(Debug, Clone)]
pub struct TidalCurrent {
    pub speed: f64,        // m/s
    pub direction: f64,    // degrees true
}

pub struct CurrentConstituent {
    pub name: String,
    pub major_axis: f64,   // m/s
    pub minor_axis: f64,   // m/s
    pub inclination: f64,  // degrees
    pub phase: f64,        // degrees
}

impl TidePredictionEngine {
    /// Predict tidal currents (simplified - real implementation needs current constituents)
    pub fn predict_current(
        &self,
        station_name: &str,
        time: DateTime<Utc>,
    ) -> Result<TidalCurrent, String> {
        // Get water level derivative as proxy for current
        let dt = Duration::minutes(30);
        let h1 = self.predict_height(station_name, time - dt)?;
        let h2 = self.predict_height(station_name, time + dt)?;
        
        // Simple relationship: current proportional to water level change
        let dh_dt = (h2 - h1) / (2.0 * dt.num_seconds() as f64 / 3600.0);
        
        // Empirical scaling (location-specific)
        let speed = (dh_dt * 0.5).abs(); // m/s
        let direction = if dh_dt > 0.0 { 45.0 } else { 225.0 }; // Flood vs ebb
        
        Ok(TidalCurrent { speed, direction })
    }
}
```

### Part 6: Applications and Visualization

```rust
use plotters::prelude::*;

/// Plot tide predictions
pub fn plot_tide_prediction(
    engine: &TidePredictionEngine,
    station_name: &str,
    start_time: DateTime<Utc>,
    days: i32,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (1200, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    // Generate predictions
    let mut predictions = Vec::new();
    let mut current_time = start_time;
    let end_time = start_time + Duration::days(days as i64);
    
    while current_time <= end_time {
        let height = engine.predict_height(station_name, current_time)?;
        predictions.push((current_time, height));
        current_time = current_time + Duration::minutes(10);
    }
    
    // Find extremes for annotations
    let extremes = engine.find_extremes(station_name, start_time, end_time)?;
    
    // Plot setup
    let min_height = predictions.iter().map(|(_, h)| h).fold(f64::INFINITY, |a, &b| a.min(b));
    let max_height = predictions.iter().map(|(_, h)| h).fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    
    let mut chart = ChartBuilder::on(&root)
        .caption(&format!("Tide Prediction - {}", station_name), ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(
            start_time..end_time,
            min_height..max_height
        )?;
    
    chart.configure_mesh()
        .x_desc("Date/Time")
        .y_desc("Water Level (m)")
        .x_label_formatter(&|x| x.format("%m/%d %H:%M").to_string())
        .draw()?;
    
    // Draw tide curve
    chart.draw_series(LineSeries::new(
        predictions.iter().map(|(t, h)| (*t, *h)),
        &BLUE,
    ))?;
    
    // Mark high and low tides
    for extreme in extremes {
        let color = if extreme.extreme_type == ExtremeType::High { &RED } else { &GREEN };
        chart.draw_series(std::iter::once(Circle::new(
            (extreme.time, extreme.height),
            5,
            color.filled(),
        )))?;
        
        // Add label
        let label = format!("{:.2}m", extreme.height);
        chart.draw_series(std::iter::once(Text::new(
            label,
            (extreme.time, extreme.height + 0.1),
            ("sans-serif", 12),
        )))?;
    }
    
    root.present()?;
    Ok(())
}

/// Create tide tables
pub fn generate_tide_table(
    engine: &TidePredictionEngine,
    station_name: &str,
    month: DateTime<Utc>,
) -> Result<String, Box<dyn std::error::Error>> {
    let start = month.with_day(1).unwrap();
    let end = if month.month() == 12 {
        start.with_year(month.year() + 1).unwrap().with_month(1).unwrap()
    } else {
        start.with_month(month.month() + 1).unwrap()
    };
    
    let extremes = engine.find_extremes(station_name, start, end)?;
    
    let mut table = format!("Tide Table - {} - {}\n", station_name, month.format("%B %Y"));
    table.push_str("Date       Time    Height  Type\n");
    table.push_str("================================\n");
    
    for extreme in extremes {
        table.push_str(&format!(
            "{} {} {:6.2}m  {}\n",
            extreme.time.format("%Y-%m-%d"),
            extreme.time.format("%H:%M"),
            extreme.height,
            match extreme.extreme_type {
                ExtremeType::High => "HIGH",
                ExtremeType::Low => "LOW ",
            }
        ));
    }
    
    Ok(table)
}
```

### Part 7: Integration with Other Coastal Processes

```rust
/// Combined water level prediction
pub struct WaterLevelPredictor {
    tide_engine: TidePredictionEngine,
    wave_predictor: Box<dyn Fn(DateTime<Utc>) -> (f64, f64)>, // Returns (Hs, setup)
    surge_predictor: Box<dyn Fn(DateTime<Utc>) -> f64>,
}

impl WaterLevelPredictor {
    /// Total water level including all components
    pub fn predict_total_level(
        &self,
        station_name: &str,
        time: DateTime<Utc>,
    ) -> Result<WaterLevelComponents, String> {
        let tide = self.tide_engine.predict_height(station_name, time)?;
        let (wave_height, wave_setup) = (self.wave_predictor)(time);
        let surge = (self.surge_predictor)(time);
        
        Ok(WaterLevelComponents {
            tide,
            surge,
            wave_setup,
            wave_height,
            total: tide + surge + wave_setup,
            extreme_level: tide + surge + wave_setup + 0.5 * wave_height, // Simplified
        })
    }
}

#[derive(Debug)]
pub struct WaterLevelComponents {
    pub tide: f64,
    pub surge: f64,
    pub wave_setup: f64,
    pub wave_height: f64,
    pub total: f64,
    pub extreme_level: f64,
}
```

### Exercises

#### Exercise 1: Tidal Classification
Implement the Form Factor to classify tidal regimes:

```rust
impl TidalStation {
    pub fn calculate_form_factor(&self) -> f64 {
        // Form Factor F = (K1 + O1) / (M2 + S2)
        // F < 0.25: Semi-diurnal
        // 0.25 < F < 1.5: Mixed semi-diurnal
        // 1.5 < F < 3.0: Mixed diurnal
        // F > 3.0: Diurnal
        // TODO: Implement
    }
}
```

#### Exercise 2: Slack Water Calculator
Find times of slack water (zero current):

```rust
pub fn find_slack_water(
    engine: &TidePredictionEngine,
    station_name: &str,
    start_time: DateTime<Utc>,
    end_time: DateTime<Utc>,
) -> Vec<DateTime<Utc>> {
    // TODO: Find when tidal current speed is minimum
    // Usually occurs near high/low tide but with lag
}
```

#### Exercise 3: Tidal Windows
Calculate operational windows for vessels:

```rust
pub struct VesselRequirements {
    pub min_depth: f64,
    pub max_current: f64,
}

pub fn find_operational_windows(
    engine: &TidePredictionEngine,
    station: &str,
    requirements: &VesselRequirements,
    date: DateTime<Utc>,
) -> Vec<(DateTime<Utc>, DateTime<Utc>)> {
    // TODO: Find periods when conditions are safe
}
```

#### Exercise 4: Real Data Integration
Fetch and analyze real tidal data:

```rust
pub async fn fetch_dfo_tide_data(
    station_id: &str,
    start_date: &str,
    end_date: &str,
) -> Result<Vec<(DateTime<Utc>, f64)>, Box<dyn std::error::Error>> {
    // TODO: Use DFO API to get water levels
    // Parse response and return time series
}
```

### Questions for Reflection

1. **Tidal Energy**: Which Canadian locations have the best tidal energy potential and why?

2. **Climate Impacts**: How will sea level rise affect tidal patterns in shallow areas like the Saint Lawrence?

3. **Ecological Timing**: How do organisms time their activities to tidal cycles?

4. **Navigation Safety**: Why are tide predictions more critical in some ports than others?

5. **Cultural Significance**: How have coastal First Nations traditionally used tidal knowledge?

### Environmental and Social Context

**Tidal Ecosystems:**
- Intertidal zones support unique biodiversity
- Tidal mixing brings nutrients
- Critical habitat for migratory birds
- Shellfish harvesting windows

**Climate Change Impacts:**
- Sea level rise changes tidal dynamics
- Altered amphidromic points
- Increased flooding during king tides
- Impacts on tidal renewable energy

**Indigenous Knowledge:**
- Traditional calendars based on tides
- Optimal harvesting times
- Navigation through tidal passages
- Oral histories of extreme events

### Integration Project: Port Safety Dashboard

Combine all Chapter 3 tutorials to create a comprehensive port safety system:

```rust
pub struct PortSafetyDashboard {
    wave_analyzer: WaveTimeSeries,
    spectral_analyzer: SpectralAnalysis,
    wave_generator: WaveGeneration,
    tide_predictor: TidePredictionEngine,
}

impl PortSafetyDashboard {
    pub fn assess_conditions(
        &self,
        port_name: &str,
        forecast_hours: i32,
    ) -> SafetyAssessment {
        // TODO: Integrate all components
        // Provide go/no-go recommendations
        // Include uncertainty estimates
    }
}
```

### Next Chapter Preview

Chapter 4 will expand our view to global wave and tidal environments, exploring:
- How BC's wave climate compares globally
- Extreme value statistics for design
- Seasonal patterns and climate variability
- Long-term trends and projections

We'll build tools to analyze decades of data and project future conditions under climate change scenarios.

### Additional Resources

- DFO Tides API: https://api-iwls.dfo-mpo.gc.ca/
- NOAA Harmonic Analysis: https://tidesandcurrents.noaa.gov/
- Indigenous Tidal Knowledge: Various coastal First Nations websites
- Tidal Power Resources: Marine Renewables Canada

Remember: Tides connect us to the cosmos - the dance of Earth, Moon, and Sun shapes our coastlines and cultures.
