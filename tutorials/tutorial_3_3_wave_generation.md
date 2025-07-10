# Tutorial 3.3: Wave Generation and Dispersion
## Ocean Waves Chapter 3.5-3.6

### Learning Objectives
By the end of this tutorial, you will:
1. Understand wind wave generation mechanisms
2. Implement wave growth models (JONSWAP, Pierson-Moskowitz)
3. Model wave dispersion and group velocity
4. Simulate wave propagation from generation to coast
5. Apply to fetch-limited conditions in Canadian waters

### Prerequisites
- Read Coastal Dynamics Ch. 3.5-3.6
- Understanding of wave dispersion relationships
- Rust development environment set up

### Introduction: From Wind to Waves

Understanding how waves are generated and propagate is crucial for:
- **Forecasting**: Predicting wave conditions from weather data
- **Hindcasting**: Reconstructing historical wave climates
- **Coastal Impact**: Understanding wave energy reaching shore
- **Navigation Safety**: Planning vessel routes

MarineLabs combines wind observations with wave physics to provide accurate forecasts, helping operators anticipate conditions hours to days in advance.

### Part 1: Wave Generation Theory

#### Key Concepts

**Energy Transfer from Wind:**
- Wind stress transfers momentum to water surface
- Phillips mechanism: Resonance with pressure fluctuations
- Miles mechanism: Shear flow instability
- Non-linear wave-wave interactions

**Fetch-Limited vs Duration-Limited:**
- Fetch: Distance over which wind blows
- Duration: Time wind has been blowing
- Fully developed sea: Energy input = dissipation

**Important Parameters:**
- Wind speed at 10m height ($U_{10}$)
- Fetch length ($F$)
- Wind duration ($t$)
- Water depth ($h$)

### Part 2: Implementing Wave Growth Models

```rust
// In src/wave_generation.rs
use std::f64::consts::PI;

// Physical constants
const G: f64 = 9.81;         // Gravitational acceleration $g$ (m/s²)
const RHO_WATER: f64 = 1025.0; // Seawater density $\rho$ (kg/m³)

/// Wave generation parameters
#[derive(Debug, Clone)]
pub struct WindConditions {
    pub speed: f64,        // m/s at 10m height
    pub direction: f64,    // degrees from North
    pub fetch: f64,        // meters
    pub duration: f64,     // seconds
    pub water_depth: f64,  // meters
}

/// Generated wave parameters
#[derive(Debug, Clone)]
pub struct GeneratedWaves {
    pub hs: f64,           // Significant wave height
    pub tp: f64,           // Peak period
    pub fp: f64,           // Peak frequency
    pub direction: f64,    // Mean wave direction
    pub is_fully_developed: bool,
}

/// Wave generation models
pub struct WaveGeneration;

impl WaveGeneration {
    /// Pierson-Moskowitz spectrum for fully developed seas
    pub fn pierson_moskowitz_hs(wind_speed: f64) -> f64 {
        // PM relationship: $H_s = 0.21 \cdot U^2 / g$
        0.21 * wind_speed.powi(2) / G
    }
    
    pub fn pierson_moskowitz_fp(wind_speed: f64) -> f64 {
        // Peak frequency: $f_p = 0.13 \cdot g / U$
        0.13 * G / wind_speed
    }
    
    /// JONSWAP fetch-limited growth
    pub fn jonswap_fetch_limited(conditions: &WindConditions) -> GeneratedWaves {
        let u = conditions.speed;
        let f = conditions.fetch;
        
        // Non-dimensional fetch
        // $\chi = g F / U^2$
        let chi = G * f / u.powi(2);
        
        // JONSWAP relationships
        let epsilon = 1.6e-3 * chi.powf(0.5);  // Non-dimensional energy: $\varepsilon = 1.6 \times 10^{-3} \chi^{0.5}$
        let nu = 3.5 * chi.powf(-0.33);        // Non-dimensional peak frequency: $\nu = 3.5 \chi^{-0.33}$
        
        // Convert to dimensional parameters
        let hs = 4.0 * (epsilon * u.powi(4) / G.powi(2)).sqrt();  // $H_s = 4\sqrt{\varepsilon U^4 / g^2}$
        let fp = nu * G / u;  // $f_p = \nu g / U$
        let tp = 1.0 / fp;
        
        // Check if fully developed
        let fully_developed_hs = Self::pierson_moskowitz_hs(u);
        let is_fully_developed = hs >= 0.95 * fully_developed_hs;
        
        GeneratedWaves {
            hs,
            tp,
            fp,
            direction: conditions.direction,
            is_fully_developed,
        }
    }
    
    /// Duration-limited growth
    pub fn duration_limited(conditions: &WindConditions) -> GeneratedWaves {
        let u = conditions.speed;
        let t = conditions.duration;
        
        // Non-dimensional time
        // $\tau = gt/U$
        let tau = G * t / u;
        
        // Duration-limited relationships
        let epsilon = 4.3e-7 * tau.powf(0.7);  // $\varepsilon = 4.3 \times 10^{-7} \tau^{0.7}$
        let nu = 0.2 * tau.powf(-0.22);  // $\nu = 0.2 \tau^{-0.22}$
        
        let hs = 4.0 * (epsilon * u.powi(4) / G.powi(2)).sqrt();
        let fp = nu * G / u;
        let tp = 1.0 / fp;
        
        GeneratedWaves {
            hs,
            tp,
            fp,
            direction: conditions.direction,
            is_fully_developed: false,
        }
    }
    
    /// Combined fetch and duration limited (takes minimum)
    pub fn fetch_duration_limited(conditions: &WindConditions) -> GeneratedWaves {
        let fetch_waves = Self::jonswap_fetch_limited(conditions);
        let duration_waves = Self::duration_limited(conditions);
        
        // Take the smaller wave height (limiting factor)
        if fetch_waves.hs < duration_waves.hs {
            fetch_waves
        } else {
            duration_waves
        }
    }
}
```

### Part 3: Wave Dispersion and Propagation

```rust
// Wave dispersion relationships
pub struct WaveDispersion;

impl WaveDispersion {
    /// Deep water dispersion relation: $\omega^2 = gk$
    pub fn deep_water_phase_speed(wavelength: f64) -> f64 {
        (G * wavelength / (2.0 * PI)).sqrt()
    }
    
    /// Finite depth dispersion relation: $\omega^2 = gk \tanh(kh)$
    pub fn finite_depth_phase_speed(wavelength: f64, depth: f64) -> f64 {
        let k = 2.0 * PI / wavelength;
        let omega_squared = G * k * (k * depth).tanh();
        let omega = omega_squared.sqrt();
        omega / k
    }
    
    /// Group velocity (energy propagation speed)
    pub fn group_velocity(wavelength: f64, depth: f64) -> f64 {
        let k = 2.0 * PI / wavelength;
        let kh = k * depth;
        // $n = \frac{1}{2}\left(1 + \frac{2kh}{\sinh(2kh)}\right)$
        let n = 0.5 * (1.0 + 2.0 * kh / kh.sinh() / kh.cosh());
        let c = Self::finite_depth_phase_speed(wavelength, depth);
        n * c  // $c_g = nc$
    }
    
    /// Check if deep water approximation is valid
    pub fn is_deep_water(wavelength: f64, depth: f64) -> bool {
        depth / wavelength > 0.5
    }
    
    /// Check if shallow water approximation is valid
    pub fn is_shallow_water(wavelength: f64, depth: f64) -> bool {
        depth / wavelength < 0.05
    }
}

/// Wave packet propagation
#[derive(Debug, Clone)]
pub struct WavePacket {
    pub position: (f64, f64),  // (x, y) coordinates
    pub energy: f64,           // Wave energy
    pub frequency: f64,        // Central frequency
    pub direction: f64,        // Propagation direction (degrees)
    pub spreading: f64,        // Directional spreading
}

impl WavePacket {
    /// Propagate wave packet over time
    pub fn propagate(&mut self, dt: f64, depth: f64) {
        // Calculate wavelength from frequency
        let wavelength = Self::wavelength_from_frequency(self.frequency, depth);
        
        // Get group velocity
        let cg = WaveDispersion::group_velocity(wavelength, depth);
        
        // Update position
        let dx = cg * dt * self.direction.to_radians().cos();
        let dy = cg * dt * self.direction.to_radians().sin();
        self.position.0 += dx;
        self.position.1 += dy;
        
        // Energy decay due to spreading (simplified)
        self.energy *= (-0.001 * dt).exp();
    }
    
    /// Iterative solution for wavelength from frequency
    fn wavelength_from_frequency(frequency: f64, depth: f64) -> f64 {
        let omega = 2.0 * PI * frequency;
        let omega_squared = omega * omega;
        
        // Initial guess (deep water)
        let mut k = omega_squared / G;
        
        // Newton-Raphson iteration
        // Solve: $\omega^2 - gk\tanh(kh) = 0$
        for _ in 0..10 {
            let f = omega_squared - G * k * (k * depth).tanh();
            let df = -G * ((k * depth).tanh() + k * depth * (1.0 - (k * depth).tanh().powi(2)));
            k -= f / df;  // $k_{n+1} = k_n - f(k_n)/f'(k_n)$
        }
        
        2.0 * PI / k
    }
}
```

### Part 4: Simulating Wave Fields

```rust
/// JONSWAP spectrum implementation
pub fn jonswap_spectrum(f: f64, hs: f64, tp: f64, gamma: f64) -> f64 {
    let fp = 1.0 / tp;
    let alpha = 0.0081; // Phillips constant
    
    // Pierson-Moskowitz part
    // $S_{PM}(f) = \alpha g^2 (2\pi)^{-4} f^{-5} \exp\left(-1.25 (f_p/f)^4\right)$
    let pm = alpha * G.powi(2) / (2.0 * PI).powi(4) / f.powi(5) 
        * (-1.25 * (fp / f).powi(4)).exp();
    
    // JONSWAP peak enhancement
    let sigma = if f <= fp { 0.07 } else { 0.09 };  // $\sigma = 0.07$ for $f \leq f_p$, $0.09$ for $f > f_p$
    let a = (-(f - fp).powi(2) / (2.0 * sigma.powi(2) * fp.powi(2))).exp();
    let enhancement = gamma.powf(a);  // $\gamma^{\exp\left(-\frac{(f-f_p)^2}{2\sigma^2 f_p^2}\right)}$
    
    // Scale to match $H_s$
    let m0_factor = hs.powi(2) / 16.0;  // $m_0 = H_s^2 / 16$
    
    pm * enhancement * m0_factor
}

/// Complete wave field simulation
pub struct WaveField {
    pub packets: Vec<WavePacket>,
    pub domain: (f64, f64, f64, f64), // (xmin, xmax, ymin, ymax)
    pub bathymetry: Box<dyn Fn(f64, f64) -> f64>, // Depth function
}

impl WaveField {
    /// Create wave field from wind conditions
    pub fn from_wind_field(
        wind_conditions: &WindConditions,
        source_location: (f64, f64),
        domain: (f64, f64, f64, f64),
        bathymetry: Box<dyn Fn(f64, f64) -> f64>,
    ) -> Self {
        let generated = WaveGeneration::fetch_duration_limited(wind_conditions);
        
        // Create spectrum of wave packets
        let mut packets = Vec::new();
        let n_frequencies = 20;
        let n_directions = 16;
        
        // Frequency range around peak
        let f_min = 0.5 * generated.fp;
        let f_max = 2.0 * generated.fp;
        
        for i in 0..n_frequencies {
            let f = f_min + (f_max - f_min) * (i as f64) / (n_frequencies as f64);
            
            // Energy from JONSWAP spectrum
            let s_f = jonswap_spectrum(f, generated.hs, generated.tp, 3.3);
            
            for j in 0..n_directions {
                let theta = generated.direction - 45.0 + 90.0 * (j as f64) / (n_directions as f64);
                
                // Directional spreading function
                let spreading = (-2.0 * ((theta - generated.direction) / 30.0).powi(2)).exp();
                
                packets.push(WavePacket {
                    position: source_location,
                    energy: s_f * spreading * (f_max - f_min) / n_frequencies as f64,
                    frequency: f,
                    direction: theta,
                    spreading: 30.0,
                });
            }
        }
        
        WaveField {
            packets,
            domain,
            bathymetry,
        }
    }
    
    /// Advance wave field in time
    pub fn advance(&mut self, dt: f64) {
        for packet in &mut self.packets {
            let depth = (self.bathymetry)(packet.position.0, packet.position.1);
            packet.propagate(dt, depth);
        }
        
        // Remove packets outside domain
        self.packets.retain(|p| {
            p.position.0 >= self.domain.0 && p.position.0 <= self.domain.1 &&
            p.position.1 >= self.domain.2 && p.position.1 <= self.domain.3
        });
    }
    
    /// Get wave statistics at a point
    pub fn get_wave_stats_at(&self, x: f64, y: f64, radius: f64) -> Option<(f64, f64)> {
        let mut total_energy = 0.0;
        let mut peak_frequency = 0.0;
        let mut max_energy = 0.0;
        
        for packet in &self.packets {
            let dx = packet.position.0 - x;
            let dy = packet.position.1 - y;
            let dist = (dx * dx + dy * dy).sqrt();
            
            if dist < radius {
                total_energy += packet.energy;
                if packet.energy > max_energy {
                    max_energy = packet.energy;
                    peak_frequency = packet.frequency;
                }
            }
        }
        
        if total_energy > 0.0 {
            let hs = 4.0 * total_energy.sqrt();
            let tp = 1.0 / peak_frequency;
            Some((hs, tp))
        } else {
            None
        }
    }
}
```

### Part 5: Application to Canadian Waters

```rust
/// Specific implementations for Canadian coastal regions
pub mod canadian_waters {
    use super::*;
    
    /// Strait of Georgia wave generation
    pub fn strait_of_georgia_waves(
        wind_speed: f64,
        wind_direction: f64,
        location: &str,
    ) -> GeneratedWaves {
        // Fetch depends on wind direction and location
        let fetch = match (location, wind_direction as i32) {
            ("Vancouver", 135..=225) => 40_000.0,  // SE winds, fetch to Gulf Islands
            ("Vancouver", 315..=45) => 100_000.0,  // NW winds, full strait fetch
            ("Nanaimo", 135..=225) => 80_000.0,    // SE winds
            _ => 50_000.0, // Default
        };
        
        let conditions = WindConditions {
            speed: wind_speed,
            direction: wind_direction,
            fetch,
            duration: 3600.0 * 6.0, // 6 hours typical
            water_depth: 200.0,     // Average depth
        };
        
        WaveGeneration::fetch_duration_limited(&conditions)
    }
    
    /// Saint Lawrence estuary waves
    pub fn st_lawrence_waves(
        wind_speed: f64,
        wind_direction: f64,
        location: &str,
    ) -> GeneratedWaves {
        // Complex fetch due to river geometry
        let (fetch, depth) = match location {
            "Quebec City" => (50_000.0, 15.0),    // Narrower, shallower
            "Tadoussac" => (100_000.0, 100.0),    // Wider, deeper
            "Sept-Iles" => (200_000.0, 50.0),     // Open to gulf
            _ => (75_000.0, 30.0),
        };
        
        let conditions = WindConditions {
            speed: wind_speed,
            direction: wind_direction,
            fetch,
            duration: 3600.0 * 12.0, // Longer duration possible
            water_depth: depth,
        };
        
        WaveGeneration::fetch_duration_limited(&conditions)
    }
}
```

### Part 6: Visualization and Analysis

```rust
use plotters::prelude::*;

pub fn visualize_wave_field(
    field: &WaveField,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Wave Field Propagation", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(
            field.domain.0..field.domain.1,
            field.domain.2..field.domain.3
        )?;
    
    chart.configure_mesh()
        .x_desc("Distance East (m)")
        .y_desc("Distance North (m)")
        .draw()?;
    
    // Draw wave packets as circles with size proportional to energy
    for packet in &field.packets {
        let size = (packet.energy * 100.0).sqrt().max(1.0);
        chart.draw_series(std::iter::once(Circle::new(
            (packet.position.0, packet.position.1),
            size,
            ShapeStyle::from(&BLUE).filled(),
        )))?;
    }
    
    root.present()?;
    Ok(())
}
```

### Exercises

#### Exercise 1: Wind Duration Calculator
Implement a function to calculate how long wind must blow to generate given wave conditions:

```rust
pub fn required_wind_duration(
    wind_speed: f64,
    target_hs: f64,
    fetch: f64,
) -> Option<f64> {
    // TODO: Inverse of duration-limited growth equations
    // Return None if target exceeds fully developed conditions
}
```

#### Exercise 2: Wave Shoaling Calculator
Model how waves change as they approach shore:

```rust
pub struct WaveShoaling {
    pub deep_water_height: f64,
    pub deep_water_period: f64,
    pub deep_water_angle: f64,
}

impl WaveShoaling {
    pub fn height_at_depth(&self, depth: f64) -> f64 {
        // TODO: Implement shoaling coefficient
        // Include refraction effects
    }
}
```

#### Exercise 3: Fetch-Limited Wave Climate
Create monthly wave statistics for a location:

```rust
pub async fn calculate_wave_climate(
    location: &str,
    year: i32,
) -> HashMap<String, WaveStatistics> {
    // TODO: Fetch historical wind data
    // Calculate wave conditions for each month
    // Return statistics
}
```

#### Exercise 4: Storm Wave Tracking
Track storm-generated waves from origin to coast:

```rust
pub struct StormTracker {
    pub storm_center: (f64, f64),
    pub max_wind_speed: f64,
    pub storm_track: Vec<(f64, f64)>,
}

impl StormTracker {
    pub fn predict_coastal_impact(
        &self,
        coastal_points: &[(f64, f64)],
    ) -> Vec<(String, f64, f64)> {
        // TODO: Generate waves along storm track
        // Propagate to coastal points
        // Return location, arrival time, and Hs
    }
}
```

### Questions for Reflection

1. **Fetch Geometry**: How does the complex coastline of BC affect wave generation compared to open ocean?

2. **Climate Change**: How might changing storm tracks affect wave climates in the Saint Lawrence?

3. **Ice Effects**: How would you modify wave generation models for partially ice-covered waters?

4. **Energy Infrastructure**: What wave conditions would be optimal for wave energy converters?

5. **Ecosystem Impacts**: How do changing wave climates affect kelp forests and eelgrass beds?

### Environmental and Social Context

**Wave Climate Changes:**
- Increasing storm intensity → higher extreme waves
- Shifting storm tracks → changed wave directions
- Reduced ice cover → longer wave generation season

**Coastal Community Impacts:**
- Erosion of culturally significant sites
- Changes to traditional harvesting areas
- Infrastructure vulnerability

**Nature-Based Solutions:**
- Kelp restoration for wave attenuation
- Living breakwaters
- Beach nourishment strategies

### Next Tutorial Preview

In Tutorial 3.4, we'll build a complete tide prediction engine, essential for:
- Port operations planning
- Intertidal habitat assessment
- Flood risk evaluation
- Traditional harvesting windows

We'll implement harmonic analysis and create predictions for Canadian coastal locations.
