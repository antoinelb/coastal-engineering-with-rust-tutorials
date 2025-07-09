mod data;
mod data_generator;
mod visualization;
mod wave_stats;

use dotenv;
use std::error::Error;
use tokio;

fn main() -> Result<(), Box<dyn Error>> {
    dotenv::dotenv()?;
    let token = dotenv::var("ONC_TOKEN").unwrap();

    let runtime = tokio::runtime::Runtime::new()?;
    runtime.block_on(data::fetch_wave_data(&token, "2024-01-01", "2024-12-31"))?;

    return Ok(());
}

// fn main() -> Result<(), Box<dyn Error>> {
//     println!("Coastal Dynamics Tutorial 3.1: Wave Statistics");
//     println!("==============================================\n");
//
//     let scenarios = [
//         ("Calm conditions", 0.5, 6.0),   // Hs=0.5m, T=6s
//         ("Moderate seas", 2.0, 8.0),     // Hs=2.0m, T=8s
//         ("Storm conditions", 4.5, 12.0), // Hs=4.5m, T=12s
//     ];
//
//     for (name, hs, period) in scenarios {
//         println!("Analyzing : {}", name);
//
//         // generate 10 minutes of data at 2 Hz
//         let data = data_generator::generate_wave_data(600.0, 2.0, hs, period);
//
//         let stats = data.calculate_statistics()?;
//         println!("  Significant wave height (Hs) : {:.2} m", stats.hs);
//         println!("  Maximum wave height (Hmax) : {:.2} m", stats.hmax);
//         println!("  Mean wave period : {:.2} s", stats.tmean);
//         println!("  Number of waves : {}", stats.wave_count);
//
//         let is_safe = data.assess_port_safety(3.0)?;
//         println!(
//             "  Sage for operations : {}",
//             if is_safe { "YES" } else { "NO" }
//         );
//
//         let p_exceed = data.exceedance_probability(3.0)?;
//         println!(
//             "  Proability of exceeding safety threshold : {:.2}%",
//             p_exceed * 100.0
//         );
//
//         let mut filename = format!(
//             "results/waves_{}.png",
//             name.to_lowercase().replace(" ", "_")
//         );
//         visualization::plot_wave_timeseries(&data, &filename)?;
//         println!("  Saved time series plot : {}", filename);
//
//         filename = format!(
//             "results/wave_heights_{}.png",
//             name.to_lowercase().replace(" ", "_")
//         );
//         visualization::plot_wave_height_distribution(
//             &data.detect_waves().unwrap(),
//             &stats,
//             &filename,
//         )?;
//         println!("  Saved wave height distribution plot : {}", filename);
//
//         println!()
//     }
//
//     println!("\nEnvironmental Context:");
//     println!("- Increased storm frequency due to climate change");
//     println!("- Wave energy impacts on coastal erosion");
//     println!("- Importance of real-time monitoring for marine safety");
//     println!("- Indigenous knowledge of seasonal wave patterns");
//
//     return Ok(());
// }
