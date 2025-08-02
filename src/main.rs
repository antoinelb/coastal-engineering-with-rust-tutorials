mod data_fetcher;
mod extreme_analysis;
mod visualization;
mod wave_climate;

use chrono::{Datelike, Duration, TimeZone, Utc};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::collections::HashMap;
use std::error::Error;
use std::f64::consts::PI;
use wave_climate::{WaveClimateAnalyzer, WaveObservation};

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {
    println!("Coastal Dynamics Tutorial 4.1: Wave Climate Analysis");
    println!("====================================================\n");

    // initialize analyzer for Tofino (west coast Vancouver Island)
    let mut analyzer = WaveClimateAnalyzer::new(
        "Tofino, BC".to_string(),
        4.0, // storm threshold: 4m
        0.5, // calm threshold: 0.5m
    );

    // generate synthetic data representing BC coast climate
    let observations = generate_bc_wave_climate(365 * 5); // 5 years
    analyzer.add_observations(observations);

    // perform analyses
    println!(
        "Analyzing {} observations...\n",
        analyzer.observations.len()
    );

    // 1. seasonal analysis
    let seasonal = analyzer.seasonal_analysis();
    println!("Seasonal Wave Climate:");
    println!(
        "Winter: mean Hs = {:.2}m, max = {:.2}m, storms = {}",
        seasonal.winter.mean_hs, seasonal.winter.max_hs, seasonal.winter.storm_count
    );
    println!(
        "Spring: mean Hs = {:.2}m, max = {:.2}m, storms = {}",
        seasonal.spring.mean_hs, seasonal.spring.max_hs, seasonal.spring.storm_count
    );
    println!(
        "Summer: mean Hs = {:.2}m, max = {:.2}m, storms = {}",
        seasonal.summer.mean_hs, seasonal.summer.max_hs, seasonal.summer.storm_count
    );
    println!(
        "Autumn: mean Hs = {:.2}m, max = {:.2}m, storms = {}",
        seasonal.autumn.mean_hs, seasonal.autumn.max_hs, seasonal.autumn.storm_count
    );

    println!();

    // 2. wave systems classification
    let (swell, sea) = analyzer.classify_wave_system();
    println!("Wave Systems:");
    println!(
        "Swell: {:.1}% of observations",
        swell.len() as f64 / analyzer.observations.len() as f64 * 100.0
    );
    println!(
        "Sea: {:.1}% of observations",
        sea.len() as f64 / analyzer.observations.len() as f64 * 100.0
    );

    println!();

    // 3. directional analysis
    let wave_rose = analyzer.wave_rose(16);
    println!(
        "Dominant Wave Direction: {:.0}°",
        wave_rose
            .iter()
            .max_by(|a, b| a.energy_fraction.partial_cmp(&b.energy_fraction).unwrap())
            .unwrap()
            .direction
    );

    println!();

    // 4. extreme value analysis
    let eva = perform_extreme_analysis(&analyzer);
    println!("Extreme Value Analysis:");
    println!("50-year return value: {:.2}m", eva.return_50);
    println!("100-year return value: {:.2}m", eva.return_100);

    // 5. generate visualizations
    visualization::plot_wave_rose(&wave_rose, "wave_rose_bc.png", "BC Coast Wave Rose")?;
    visualization::plot_seasonal_comparison(&analyzer, "seasonal_bc.png")?;

    Ok(())
}

// generate synthetic wave dat representative of BC coast
fn generate_bc_wave_climate(days: usize) -> Vec<WaveObservation> {
    let mut rng = rand::rng();
    let mut observations = Vec::new();
    let start_date = Utc.with_ymd_and_hms(2019, 1, 1, 0, 0, 0).unwrap();

    for day in 0..days {
        let date = start_date + Duration::days(day as i64);
        let day_of_year = date.ordinal();

        // seasonal variation (stronger in winter)
        let seasonal_factor = 1.0 + 0.8 + (2.0 * PI * (day_of_year as f64 - 80.0) / 365.0).cos();

        // generate 8 observations per day (3-hour intervals)
        for hour in (0..24).step_by(3) {
            let timestamp = date + Duration::hours(hour as i64);

            // base conditions with seasonal variation
            let base_hs = 1.5 * seasonal_factor;
            let hs_dist = Normal::new(base_hs, 0.5 * seasonal_factor).unwrap();
            let hs = hs_dist.sample(&mut rng).max(0.1);

            // period correlates with height
            let tp = 8.0 + 2.0 * hs.sqrt() + rng.random::<f64>() * 2.0;
            let tm = tp * (0.7 + rng.random::<f64>() * 0.2);

            // direction: mainly from W-NW (270-315°) with some spread
            let dp = 292.5 + Normal::new(0.0, 20.0).unwrap().sample(&mut rng);
            let dm = dp + Normal::new(0.0, 5.0).unwrap().sample(&mut rng);

            // wind correlation
            let wind_speed = Some(5.0 + 3.0 * hs + rng.random::<f64>() * 5.0);
            let wind_dir = Some(dp + Normal::new(-10.0, 15.0).unwrap().sample(&mut rng));

            observations.push(WaveObservation {
                timestamp,
                hs,
                tp,
                tm,
                dp: dp % 360.0,
                dm: dm % 360.0,
                spread: 25.0 + rng.random::<f64>() * 10.0,
                wind_speed,
                wind_dir,
            });
        }
    }

    observations
}

struct ExtremeAnalysisResults {
    return_50: f64,
    return_100: f64,
}

fn perform_extreme_analysis(analyzer: &WaveClimateAnalyzer) -> ExtremeAnalysisResults {
    let mut annual_maxima: HashMap<i32, f64> = HashMap::new();

    for obs in &analyzer.observations {
        let year = obs.timestamp.year();
        let entry = annual_maxima.entry(year).or_insert(0.0);
        if obs.hs > *entry {
            *entry = obs.hs;
        }
    }

    let maxima: Vec<f64> = annual_maxima.values().cloned().collect();

    // simple Gumbel fit (assuming GEV shape parameter ≈ 0)
    let mean = maxima.iter().sum::<f64>() / maxima.len() as f64;
    let variance = maxima.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / maxima.len() as f64;
    let std = variance.sqrt();

    let scale = std * (6.0_f64).sqrt() / PI;
    let location = mean - 0.5772 * scale; // Euler's constant

    // return values using Gumbel distribution
    let return_50 = location - scale * (-(50.0_f64.ln())).ln();
    let return_100 = location - scale * (-(100.0_f64.ln())).ln();

    ExtremeAnalysisResults {
        return_50,
        return_100,
    }
}
