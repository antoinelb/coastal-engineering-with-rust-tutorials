use crate::wave_generation::{GeneratedWaves, WaveGeneration, WindConditions};

pub fn strait_of_georgia_waves(
    wind_speed: f64,
    wind_direction: f64,
    location: &str,
) -> GeneratedWaves {
    // fetch depends on wind direction and location
    let fetch = match (location, wind_direction as i32) {
        ("Vancouver", 135..=225) => 40_000.0, // SE winds, fetch to Gulf Islands
        ("Vancouver", 45..=315) => 100_000.0, // NW winds, full strait fetch
        ("Nanaimo", 135..=225) => 80_000.0,   // SE winds
        _ => 50_000.0,
    };

    let conditions = WindConditions {
        speed: wind_speed,
        direction: wind_direction,
        fetch,
        duration: 3600.0 * 6.0, // 6 hours typical
        water_depth: 200.0,     // average depth
    };

    return WaveGeneration::fetch_duration_limited(&conditions);
}

pub fn st_lawrence_waves(wind_speed: f64, wind_direction: f64, location: &str) -> GeneratedWaves {
    // complex fetch due to river geometry
    let (fetch, depth) = match location {
        "Quebec City" => (50_000.0, 15.0), // narrower, shallower
        "Tadoussac" => (100_000.0, 100.0), // wider, deeper
        "Sept-Iles" => (200_000.0, 50.0),  // open to gulf
        _ => (75_000.0, 30.0),
    };

    let conditions = WindConditions {
        speed: wind_speed,
        direction: wind_direction,
        fetch,
        duration: 3600.0 * 6.0, // longer duration possible
        water_depth: depth,
    };

    return WaveGeneration::fetch_duration_limited(&conditions);
}
