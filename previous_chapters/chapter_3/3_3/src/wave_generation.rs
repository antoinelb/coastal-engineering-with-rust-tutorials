use std::f64::consts::PI;

const G: f64 = 9.81;
const RHO_WATER: f64 = 1025.0;

#[derive(Debug, Clone)]
pub struct WindConditions {
    pub speed: f64,       // m/s at 10m height
    pub direction: f64,   // degrees from north
    pub fetch: f64,       // meters
    pub duration: f64,    // seconds
    pub water_depth: f64, // meters
}

pub struct GeneratedWaves {
    pub hs: f64,        // significant wave height
    pub tp: f64,        // peak period
    pub fp: f64,        // peak frequency
    pub direction: f64, // mean wave direction
    pub is_fully_developed: bool,
}

pub struct WaveGeneration;

impl WaveGeneration {
    // Pierson-Moskowitz spectrum for fully developed seas
    pub fn pierson_moskowitz_hs(wind_speed: f64) -> f64 {
        return 0.21 * wind_speed.powi(2) / G;
    }

    pub fn pierson_moskowitz_fp(wind_speed: f64) -> f64 {
        return 0.13 * G / wind_speed;
    }

    // JONSWAP fetch-limited growth
    pub fn jonswap_fetch_limited(conditions: &WindConditions) -> GeneratedWaves {
        let u = conditions.speed;
        let f = conditions.fetch;

        // non-dimensional fetch
        let chi = G * f / u.powi(2);

        // JONSWAP relationships
        let epsilon = 1.6e-3 * chi.powf(0.5); // non-dimensional energy
        let nu = 3.5 * chi.powf(-0.33); // non-dimensional peak frequency

        // convert to dimensional parameters
        let hs = 4.0 * (epsilon * u.powi(4) / G.powi(2)).sqrt();
        let fp = nu * G / u;
        let tp = 1.0 / fp;

        // check if fully developed
        let fully_developed_hs = Self::pierson_moskowitz_hs(u);
        let is_fully_developed = hs >= 0.95 * fully_developed_hs;

        return GeneratedWaves {
            hs,
            tp,
            fp,
            direction: conditions.direction,
            is_fully_developed,
        };
    }

    // duration-limited growth
    pub fn duration_limited(conditions: &WindConditions) -> GeneratedWaves {
        let u = conditions.speed;
        let t = conditions.duration;

        // non-dimensional time
        let tau = G * t / u;

        // duration-limited relationships
        let epsilon = 4.3e-7 * tau.powf(0.7);
        let nu = 0.2 * tau.powf(-0.22);
        let hs = 4.0 * (epsilon * u.powi(4) / G.powi(2)).sqrt();
        let fp = nu * G / u;
        let tp = 1.0 / fp;

        return GeneratedWaves {
            hs,
            tp,
            fp,
            direction: conditions.direction,
            is_fully_developed: false,
        };
    }

    // combined fetch and duration limited (minimum of the two)
    pub fn fetch_duration_limited(conditions: &WindConditions) -> GeneratedWaves {
        let fetch_waves = Self::jonswap_fetch_limited(conditions);
        let duration_waves = Self::duration_limited(conditions);

        if fetch_waves.hs < duration_waves.hs {
            return fetch_waves;
        } else {
            return duration_waves;
        }
    }
}

pub struct WaveDispersion;

impl WaveDispersion {
    // deep wawter dispersion relation
    pub fn deep_Water_phase_speed(wavelength: f64) -> f64 {
        return (G * wavelength / (2.0 * PI)).sqrt();
    }

    // finite depth dispersion relation
    pub fn finite_depth_phase_speed(wavelength: f64, depth: f64) -> f64 {
        let k = 2.0 * PI / wavelength;
        let omega_squared = G * k * (k * depth).tanh();
        return omega_squared.sqrt() / k;
    }

    // group velocity (energy propagation speed)
    pub fn group_velocity(wavelength: f64, depth: f64) -> f64 {
        let k = 2.0 * PI / wavelength;
        let kh = k * depth;
        let n = 0.5 * (1.0 + 2.0 * kh / kh.sinh() / kh.cos());
        let c = Self::finite_depth_phase_speed(wavelength, depth);
        return n * c;
    }

    // check if deep water approximation is valid
    pub fn is_deep_water(wavelength: f64, depth: f64) -> bool {
        return depth / wavelength > 0.5;
    }

    // check if shallow water approximation is valid
    pub fn is_shallow_water(wavelength: f64, depth: f64) -> bool {
        return depth / wavelength < 0.05;
    }
}

// wave packet propagation
#[derive(Debug, Clone)]
pub struct WavePacket {
    pub position: (f64, f64), // (x, y) coordinates
    pub energy: f64,          // wave energy
    pub frequency: f64,       // central frequency
    pub direction: f64,       // propagation direction (degrees)
    pub spreading: f64,       // directional spreading
}

impl WavePacket {
    // propagate wave packet over time
    pub fn propagate(&mut self, dt: f64, depth: f64) {
        // calculate wavelength from frequency
        let wavelength = Self::wavelength_from_frequency(self.frequency, depth);

        // get group velocity
        let cg = WaveDispersion::group_velocity(wavelength, depth);

        // update position
        let dx = cg * dt * self.direction.to_radians().cos();
        let dy = cg * dt * self.direction.to_radians().sin();
        self.position.0 += dx;
        self.position.1 += dy;

        // energy decay due to spreading (simplified)
        self.energy *= (-0.001 * dt).exp();
    }

    // iterative solution for wavelength from frequency
    fn wavelength_from_frequency(frequency: f64, depth: f64) -> f64 {
        let omega = 2.0 * PI * frequency;
        let omega_squared = omega * omega;

        // deep water as initial guess
        let mut k = omega_squared / G;

        // newton-raphson iterations
        for _ in 0..10 {
            let f = omega_squared - G * k * (k * depth).tanh();
            let df = -G * ((k * depth).tanh() + k * depth * (1.0 - (k * depth).tanh().powi(2)));
            k -= f / df;
        }

        return 2.0 * PI / k;
    }
}

// JONSWAP spectrum implementation
pub fn jonswap_sepectrum(f: f64, hs: f64, tp: f64, gamma: f64) -> f64 {
    let fp = 1.0 / tp;
    let alpha = 0.0081; // Phillips constant

    // Pierson-Moskowitz part
    let pm = alpha * G.powi(2) / (2.0 * PI).powi(4) / f.powi(5) * (-1.25 * (fp / f).powi(4)).exp();

    // JONSWAP peak enhancement
    let sigma: f64 = if f <= fp { 0.07 } else { 0.09 };
    let a = (-(f - fp).powi(2) / (2.0 * sigma.powi(2) * fp.powi(2))).exp();
    let enhancement = gamma.powf(a);

    // scale to match Hs
    let m0_factor = hs.powi(2) / 16.0;

    return pm * enhancement * m0_factor;
}

pub struct WaveField {
    pub packets: Vec<WavePacket>,
    pub domain: (f64, f64, f64, f64), // (xmin, xmax, ymin, ymax)
    pub bathymetry: Box<dyn Fn(f64, f64) -> f64>, // depth function
}

impl WaveField {
    // create wave field from wind conditions
    pub fn from_wind_field(
        wind_conditions: &WindConditions,
        source_location: (f64, f64),
        domain: (f64, f64, f64, f64),
        bathymetry: Box<dyn Fn(f64, f64) -> f64>,
    ) -> Self {
        let generated = WaveGeneration::fetch_duration_limited(wind_conditions);

        // create spectrum of wind packets
        let mut packets = Vec::new();
        let n_frequencies = 20;
        let n_directions = 16;

        // frequency range around peak
        let f_min = 0.5 * generated.fp;
        let f_max = 2.0 * generated.fp;

        for i in 0..n_frequencies {
            let f = f_min + (f_max - f_min) * (i as f64) / (n_frequencies as f64);

            // energy from jonwap spectrum
            let s_f = jonswap_sepectrum(f, generated.hs, generated.tp, 3.3);

            for j in 0..n_directions {
                let theta = generated.direction - 45.0 + 90.0 * (j as f64) / (n_directions as f64);

                // directional spreading function
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

        return WaveField {
            packets,
            domain,
            bathymetry,
        };
    }

    // advance wave field in time
    pub fn advance(&mut self, dt: f64) {
        for packet in &mut self.packets {
            let depth = (self.bathymetry)(packet.position.0, packet.position.1);
            packet.propagate(dt, depth);
        }

        // remove packets outside domain
        self.packets.retain(|p| {
            p.position.0 >= self.domain.0
                && p.position.0 <= self.domain.1
                && p.position.1 >= self.domain.2
                && p.position.1 <= self.domain.3
        });
    }

    // get wave statistics at a point
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
            return Some((hs, tp));
        } else {
            return None;
        }
    }
}
