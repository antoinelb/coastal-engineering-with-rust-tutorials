# Tutorial 3.1: Wave Statistics and Rust Fundamentals
## Ocean Waves Chapter 3.1-3.4

### Learning Objectives
By the end of this tutorial, you will:
1. Understand wave height distributions and statistical descriptors
2. Implement basic wave analysis functions in Rust
3. Learn Rust fundamentals: ownership, borrowing, and error handling
4. Process time series data using Rust vectors and iterators
5. Visualize wave statistics relevant to coastal safety

### Prerequisites
- Read Coastal Dynamics Ch. 3.1-3.4
- Basic understanding of wave mechanics
- Rust development environment set up

### Introduction: Waves and Marine Safety

MarineLabs deploys sensor networks to measure waves in real-time, providing critical data for port operations and coastal safety. Understanding wave statistics is fundamental to:
- Predicting extreme events (like the 17.6m rogue wave near Ucluelet)
- Optimizing vessel operations during rough conditions
- Designing coastal infrastructure

In this tutorial, we'll build the foundation for a wave analysis system similar to what powers MarineLabs' CoastAware platform.

### Part 1: Understanding Wave Measurements

#### Key Concepts Review

**Wave Parameters:**
- Wave height ($H$): Vertical distance from trough to crest
- Wave period ($T$): Time between successive crests
- Significant wave height ($H_s$): Average of highest 1/3 of waves
- Maximum wave height ($H_{max}$): Highest individual wave

**Why These Matter:**
- Port operations suspend when $H_s$ exceeds safety thresholds
- Vessel scheduling depends on wave period (longer periods = safer conditions)
- Extreme waves pose risks to infrastructure and vessels

#### Questions to Consider
1. Why is significant wave height ($H_s$) more useful than mean wave height for marine operations?
2. How might climate change affect wave statistics along the BC coast?
3. What environmental factors influence wave generation in the Strait of Georgia?

### Part 2: Rust Fundamentals Through Wave Analysis

Let's start building our wave analysis toolkit. Create a new file `src/wave_stats.rs`:

```rust
// wave_stats.rs
use std::error::Error;
use std::fmt;

/// Custom error type for wave analysis
#[derive(Debug)]
pub enum WaveError {
    InsufficientData,
    InvalidMeasurement(String),
}

impl fmt::Display for WaveError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            WaveError::InsufficientData => write!(f, "Insufficient wave data"),
            WaveError::InvalidMeasurement(msg) => write!(f, "Invalid measurement: {}", msg),
        }
    }
}

impl Error for WaveError {}

/// Represents a single wave measurement
#[derive(Debug, Clone)]
pub struct WaveMeasurement {
    pub time: f64,      // seconds since start
    pub elevation: f64, // meters
}

/// Container for wave time series data
pub struct WaveTimeSeries {
    measurements: Vec<WaveMeasurement>,
    sampling_rate: f64, // Hz
}

impl WaveTimeSeries {
    /// Create a new wave time series
    pub fn new(sampling_rate: f64) -> Self {
        WaveTimeSeries {
            measurements: Vec::new(),
            sampling_rate,
        }
    }

    /// Add a measurement (demonstrates ownership)
    pub fn add_measurement(&mut self, measurement: WaveMeasurement) {
        self.measurements.push(measurement);
    }

    /// Get number of measurements (demonstrates borrowing)
    pub fn len(&self) -> usize {
        self.measurements.len()
    }

    /// Calculate mean water level
    pub fn mean_level(&self) -> Result<f64, WaveError> {
        if self.measurements.is_empty() {
            return Err(WaveError::InsufficientData);
        }

        let sum: f64 = self.measurements
            .iter()
            .map(|m| m.elevation)
            .sum();

        Ok(sum / self.measurements.len() as f64)
    }
}
```

**Rust Concepts Demonstrated:**
1. **Ownership**: `add_measurement` takes ownership of the measurement
2. **Borrowing**: `&self` and `&mut self` borrow the struct
3. **Error Handling**: Using `Result<T, E>` for fallible operations
4. **Iterators**: Using `iter()` and `map()` for functional programming

### Part 3: Implementing Zero-Crossing Analysis

Zero-crossing analysis is fundamental to wave statistics. Add this to `wave_stats.rs`:

```rust
/// Represents an individual wave
#[derive(Debug, Clone)]
pub struct Wave {
    pub height: f64,    // meters
    pub period: f64,    // seconds
    pub crest: f64,     // maximum elevation
    pub trough: f64,    // minimum elevation
}

impl WaveTimeSeries {
    /// Detect individual waves using zero-crossing method
    pub fn detect_waves(&self) -> Result<Vec<Wave>, WaveError> {
        if self.measurements.len() < 3 {
            return Err(WaveError::InsufficientData);
        }

        // Remove mean level
        let mean = self.mean_level()?;
        let detrended: Vec<f64> = self.measurements
            .iter()
            .map(|m| m.elevation - mean)
            .collect();

        let mut waves = Vec::new();
        let mut i = 0;

        // Find zero crossings
        while i < detrended.len() - 1 {
            // Find upward zero crossing
            if detrended[i] <= 0.0 && detrended[i + 1] > 0.0 {
                let start_idx = i;
                
                // Find next downward zero crossing
                while i < detrended.len() - 1 && detrended[i + 1] >= 0.0 {
                    i += 1;
                }
                
                // Find next upward zero crossing (end of wave)
                while i < detrended.len() - 1 && detrended[i + 1] < 0.0 {
                    i += 1;
                }
                
                if i < detrended.len() - 1 {
                    let end_idx = i;
                    
                    // Extract wave properties
                    let wave_segment = &self.measurements[start_idx..=end_idx];
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
                    
                    waves.push(Wave {
                        height,
                        period,
                        crest,
                        trough,
                    });
                }
            }
            i += 1;
        }

        if waves.is_empty() {
            return Err(WaveError::InsufficientData);
        }

        Ok(waves)
    }
}
```

### Part 4: Calculating Wave Statistics

Now let's implement the key statistical measures:

```rust
/// Wave statistics summary
#[derive(Debug)]
pub struct WaveStatistics {
    pub hs: f64,           // Significant wave height
    pub hmax: f64,         // Maximum wave height
    pub hmean: f64,        // Mean wave height
    pub tmean: f64,        // Mean wave period
    pub wave_count: usize, // Number of waves
}

impl WaveTimeSeries {
    /// Calculate comprehensive wave statistics
    pub fn calculate_statistics(&self) -> Result<WaveStatistics, WaveError> {
        let waves = self.detect_waves()?;
        
        if waves.is_empty() {
            return Err(WaveError::InsufficientData);
        }

        // Sort waves by height (descending)
        let mut sorted_waves = waves.clone();
        sorted_waves.sort_by(|a, b| b.height.partial_cmp(&a.height).unwrap());

        // Calculate $H_s$ (average of highest 1/3)
        let n_third = (sorted_waves.len() as f64 / 3.0).ceil() as usize;
        let hs = sorted_waves[..n_third]
            .iter()
            .map(|w| w.height)
            .sum::<f64>() / n_third as f64;

        // Other statistics
        let hmax = sorted_waves[0].height;  // $H_{max}$
        let hmean = waves.iter().map(|w| w.height).sum::<f64>() / waves.len() as f64;  // $H_{mean}$
        let tmean = waves.iter().map(|w| w.period).sum::<f64>() / waves.len() as f64;  // $T_{mean}$

        Ok(WaveStatistics {
            hs,
            hmax,
            hmean,
            tmean,
            wave_count: waves.len(),
        })
    }

    /// Check if conditions are safe for port operations
    pub fn assess_safety(&self, hs_threshold: f64) -> Result<bool, WaveError> {
        let stats = self.calculate_statistics()?;
        Ok(stats.hs < hs_threshold)
    }
}
```

### Part 5: Working with Real Data

Let's create a data generator that simulates realistic wave conditions:

```rust
// In src/data_generator.rs
use rand::{thread_rng, Rng};
use rand_distr::{Normal, Distribution};
use crate::wave_stats::{WaveMeasurement, WaveTimeSeries};

/// Generate synthetic wave data similar to BC coast conditions
pub fn generate_wave_data(
    duration: f64,      // seconds
    sampling_rate: f64, // Hz
    hs: f64,           // target significant wave height
    peak_period: f64,  // dominant wave period
) -> WaveTimeSeries {
    let mut rng = thread_rng();
    let mut time_series = WaveTimeSeries::new(sampling_rate);
    
    let dt = 1.0 / sampling_rate;
    let n_samples = (duration * sampling_rate) as usize;
    
    // Generate using superposition of sinusoids (simplified)
    for i in 0..n_samples {
        let t = i as f64 * dt;
        let mut elevation = 0.0;
        
        // Add primary wave component
        let primary_amplitude = hs * 0.7;
        let primary_frequency = 1.0 / peak_period;  // $f = 1/T$
        elevation += primary_amplitude * (2.0 * std::f64::consts::PI * primary_frequency * t).sin();  // $A\sin(2\pi f t)$
        
        // Add secondary components
        for j in 1..5 {
            let amplitude = primary_amplitude / (j as f64 * 1.5);
            let frequency = primary_frequency * (0.8 + 0.1 * j as f64);
            let phase = rng.gen::<f64>() * 2.0 * std::f64::consts::PI;
            elevation += amplitude * (2.0 * std::f64::consts::PI * frequency * t + phase).sin();
        }
        
        // Add noise
        let noise_dist = Normal::new(0.0, hs * 0.1).unwrap();
        elevation += noise_dist.sample(&mut rng);
        
        time_series.add_measurement(WaveMeasurement { time: t, elevation });
    }
    
    time_series
}
```

### Part 6: Visualization

Create visualization functions using plotters:

```rust
// In src/visualization.rs
use plotters::prelude::*;
use crate::wave_stats::{WaveTimeSeries, WaveStatistics};

pub fn plot_wave_timeseries(
    data: &WaveTimeSeries,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let measurements: Vec<(f64, f64)> = data.measurements
        .iter()
        .map(|m| (m.time, m.elevation))
        .collect();
    
    let max_time = measurements.last().map(|(t, _)| *t).unwrap_or(0.0);
    let (min_elev, max_elev) = measurements
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), (_, e)| {
            (min.min(*e), max.max(*e))
        });
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Wave Elevation Time Series", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..max_time, min_elev..max_elev)?;
    
    chart.configure_mesh()
        .x_desc("Time (s)")
        .y_desc("Elevation (m)")
        .draw()?;
    
    chart.draw_series(LineSeries::new(measurements, &BLUE))?
        .label("Surface elevation")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &BLUE));
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
    
    root.present()?;
    Ok(())
}

pub fn plot_wave_height_distribution(
    waves: &[Wave],
    stats: &WaveStatistics,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    // Create histogram
    let max_height = waves.iter().map(|w| w.height).fold(0.0, f64::max);
    let bin_width = max_height / 20.0;
    let mut histogram = vec![0; 20];
    
    for wave in waves {
        let bin = ((wave.height / bin_width) as usize).min(19);
        histogram[bin] += 1;
    }
    
    let max_count = *histogram.iter().max().unwrap() as f64;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Wave Height Distribution", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..max_height, 0.0..max_count)?;
    
    chart.configure_mesh()
        .x_desc("Wave Height (m)")
        .y_desc("Count")
        .draw()?;
    
    // Draw histogram bars
    chart.draw_series(
        histogram.iter().enumerate().map(|(i, &count)| {
            Rectangle::new([
                (i as f64 * bin_width, 0.0),
                ((i + 1) as f64 * bin_width, count as f64)
            ], BLUE.filled())
        })
    )?;
    
    // Add vertical lines for statistics
    chart.draw_series(vec![
        PathElement::new(vec![(stats.hs, 0.0), (stats.hs, max_count)], &RED),
        PathElement::new(vec![(stats.hmax, 0.0), (stats.hmax, max_count)], &GREEN),
    ])?;
    
    // Add legend
    chart.draw_series(std::iter::once(PathElement::new(vec![(0.0, 0.0)], &RED)))?
        .label(format!("Hs = {:.2} m", stats.hs))
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &RED));
    
    chart.draw_series(std::iter::once(PathElement::new(vec![(0.0, 0.0)], &GREEN)))?
        .label(format!("Hmax = {:.2} m", stats.hmax))
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &GREEN));
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
    
    root.present()?;
    Ok(())
}
```

### Part 7: Putting It All Together

Create the main application in `src/main.rs`:

```rust
mod wave_stats;
mod data_generator;
mod visualization;

use wave_stats::WaveTimeSeries;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    println!("Coastal Dynamics Tutorial 3.1: Wave Statistics");
    println!("=============================================\n");

    // Simulate different sea states
    let scenarios = vec![
        ("Calm conditions", 0.5, 6.0),    // $H_s = 0.5\text{m}$, $T = 6\text{s}$
        ("Moderate seas", 2.0, 8.0),       // $H_s = 2.0\text{m}$, $T = 8\text{s}$
        ("Storm conditions", 4.5, 12.0),   // $H_s = 4.5\text{m}$, $T = 12\text{s}$
    ];

    for (name, hs, period) in scenarios {
        println!("Analyzing: {}", name);
        
        // Generate 10 minutes of data at 2 Hz
        let data = data_generator::generate_wave_data(600.0, 2.0, hs, period);
        
        // Calculate statistics
        let stats = data.calculate_statistics()?;
        println!("  Significant wave height (Hs): {:.2} m", stats.hs);
        println!("  Maximum wave height (Hmax): {:.2} m", stats.hmax);
        println!("  Mean wave period: {:.2} s", stats.tmean);
        println!("  Number of waves: {}", stats.wave_count);
        
        // Safety assessment (typical BC ferry threshold)
        let safe = data.assess_safety(3.0)?;
        println!("  Safe for operations: {}", if safe { "YES" } else { "NO" });
        
        // Visualize
        let filename = format!("waves_{}.png", name.replace(" ", "_"));
        visualization::plot_wave_timeseries(&data, &filename)?;
        println!("  Saved time series plot: {}", filename);
        
        println!();
    }

    // Environmental considerations
    println!("\nEnvironmental Context:");
    println!("- Increased storm frequency due to climate change");
    println!("- Wave energy impacts on coastal erosion");
    println!("- Importance of real-time monitoring for marine safety");
    println!("- Indigenous knowledge of seasonal wave patterns");

    Ok(())
}
```

### Exercises

#### Exercise 1: Rayleigh Distribution
The Rayleigh distribution is commonly used to model wave heights. Implement:

```rust
impl WaveStatistics {
    /// Calculate Rayleigh distribution parameters
    pub fn rayleigh_parameters(&self) -> (f64, f64) {
        // TODO: Calculate scale parameter from $H_{rms}$
        // $H_{rms} = \sqrt{\overline{H^2}}$
        // $\text{scale} = H_{rms} / \sqrt{\pi/2}$
    }
    
    /// Probability of exceeding given wave height
    pub fn exceedance_probability(&self, height: f64) -> f64 {
        // TODO: Use Rayleigh CDF
        // $P(H > h) = \exp\left(-\frac{h^2}{2\sigma^2}\right)$
    }
}
```

#### Exercise 2: Wave Climate Analysis
Create a function to analyze seasonal variations:

```rust
pub struct SeasonalWaveClimate {
    pub winter_hs: f64,
    pub summer_hs: f64,
    pub storm_frequency: f64,
}

impl SeasonalWaveClimate {
    pub fn from_annual_data(/* TODO: parameters */) -> Self {
        // Analyze wave data by season
        // Consider BC coast patterns
    }
}
```

#### Exercise 3: Real Data Integration
Implement a function to fetch real wave data from Ocean Networks Canada:

```rust
async fn fetch_wave_data(station: &str, start_date: &str, end_date: &str) 
    -> Result<WaveTimeSeries, Box<dyn Error>> {
    // TODO: Use reqwest to fetch data
    // Parse JSON response
    // Convert to WaveTimeSeries
}
```

#### Exercise 4: Port Safety Dashboard
Create a simple terminal-based dashboard showing:
- Current sea state
- Safety status for different vessel types
- Forecast trends
- Historical extremes

### Questions for Reflection

1. **Climate Impact**: How might changing wave climates affect coastal communities along the Saint-Lawrence?

2. **Engineering Design**: If designing a harbour entrance, what wave statistics would be most critical?

3. **Environmental Justice**: How do wave-induced coastal changes disproportionately affect certain communities?

4. **Data Quality**: What challenges might sensor networks face in harsh marine environments?

5. **Rust Design**: Why is Rust's ownership model particularly suitable for real-time data processing?

### Additional Resources

- Ocean Networks Canada API documentation
- DFO wave measurement standards
- Rust async programming guide
- MarineLabs technical papers on wave measurement

### Next Tutorial Preview

In Tutorial 3.2, we'll implement spectral analysis to understand wave energy distribution across frequencies - a key tool for predicting vessel response and coastal impacts.

Remember: The ocean is a complex system. Our models are simplifications, but they provide valuable insights for safety and sustainability.
