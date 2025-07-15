mod tide_prediction;

use chrono::{DateTime, Duration, Utc};
use tide_prediction::{TidalConstituent, TidalStation, ExtremeType};

fn main() {
    // Example usage of TidalStation public methods
    
    // Create a sample tidal station with major constituents
    let station = TidalStation {
        name: "Halifax Harbour".to_string(),
        lat: 44.6488,
        lon: -63.5752,
        timezone: -4,
        datum_offset: 0.0,
        constituents: vec![
            TidalConstituent {
                name: "M2".to_string(),
                frequency: 1.9323,
                amplitude: 0.98,
                phase: 123.4,
                speed: 28.9841,
            },
            TidalConstituent {
                name: "S2".to_string(),
                frequency: 2.0,
                amplitude: 0.31,
                phase: 156.7,
                speed: 30.0,
            },
            TidalConstituent {
                name: "O1".to_string(),
                frequency: 0.9295,
                amplitude: 0.15,
                phase: 78.2,
                speed: 13.9430,
            },
        ],
    };
    
    // Example 1: predict_height() - Get tide height at a specific time
    let current_time = Utc::now();
    match station.predict_height(current_time) {
        Ok(height) => println!("Current tide height: {:.2} m", height),
        Err(e) => println!("Error predicting height: {}", e),
    }
    
    // Example 2: find_extremes() - Find high and low tides over a 24-hour period
    let start_time = current_time;
    let end_time = start_time + Duration::hours(24);
    
    match station.find_extremes(start_time, end_time) {
        Ok(extremes) => {
            println!("\nTide extremes for next 24 hours:");
            for extreme in extremes {
                let tide_type = match extreme.extreme_type {
                    ExtremeType::High => "High",
                    ExtremeType::Low => "Low",
                };
                println!("{} tide at {}: {:.2} m", 
                    tide_type, 
                    extreme.time.format("%Y-%m-%d %H:%M UTC"),
                    extreme.height
                );
            }
        },
        Err(e) => println!("Error finding extremes: {}", e),
    }
}
