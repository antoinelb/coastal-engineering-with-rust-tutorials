use crate::waves::{WaveMeasurement, WaveTimeSeries};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::f64::consts::PI;

pub fn generate_wave_data(
    duration: f64,
    sampling_rate: f64,
    hs: f64,
    peak_period: f64,
) -> WaveTimeSeries {
    let mut rng = rand::rng();
    let mut time_series = WaveTimeSeries::new(sampling_rate);

    let dt = 1.0 / sampling_rate;
    let n_samples = (duration * sampling_rate) as usize;

    // generate using superposition of sinusoids
    for i in 0..n_samples {
        let t = i as f64 * dt;
        let mut elevation = 0.0;

        // primary wave component
        let primary_amplitude = hs * 0.7;
        let primary_frequency = 1.0 / peak_period;
        elevation += primary_amplitude * (2.0 * PI * primary_frequency * t).sin();

        // secondary wave components
        for j in 1..5 {
            let amplitude = primary_amplitude / (j as f64 * 1.5);
            let frequency = primary_frequency * (0.8 + 0.1 * j as f64);
            let phase = rng.random::<f64>() * 2.0 * PI;
            elevation += amplitude * (2.0 * PI * frequency * t + phase).sin();
        }

        // add noise
        let noise_dist = Normal::new(0.0, hs * 1.0).unwrap();
        elevation += noise_dist.sample(&mut rng);

        time_series.add_measurement(WaveMeasurement { time: t, elevation });
    }

    return time_series;
}
