mod canadian_waters;
mod visualization;
mod wave_generation;

use canadian_waters::{strait_of_georgia_waves, st_lawrence_waves};
use wave_generation::{
    WaveGeneration, WaveDispersion, WavePacket, WaveField, WindConditions, 
    jonswap_sepectrum
};
use visualization::visualize_wave_field;

fn main() {
    println!("=== Coastal Dynamics Simulation Toolkit ===\n");
    
    // Demonstrate Canadian Waters module
    println!("1. CANADIAN WATERS MODULE");
    println!("========================");
    
    // Strait of Georgia scenarios
    println!("\n1.1 Strait of Georgia Waves:");
    let vancouver_se = strait_of_georgia_waves(15.0, 180.0, "Vancouver");
    println!("Vancouver (15 m/s SE wind): Hs = {:.2} m, Tp = {:.2} s", 
             vancouver_se.hs, vancouver_se.tp);
    
    let vancouver_nw = strait_of_georgia_waves(20.0, 315.0, "Vancouver");
    println!("Vancouver (20 m/s NW wind): Hs = {:.2} m, Tp = {:.2} s", 
             vancouver_nw.hs, vancouver_nw.tp);
    
    let nanaimo_se = strait_of_georgia_waves(12.0, 180.0, "Nanaimo");
    println!("Nanaimo (12 m/s SE wind): Hs = {:.2} m, Tp = {:.2} s", 
             nanaimo_se.hs, nanaimo_se.tp);
    
    // St. Lawrence scenarios
    println!("\n1.2 St. Lawrence River Waves:");
    let quebec_city = st_lawrence_waves(18.0, 90.0, "Quebec City");
    println!("Quebec City (18 m/s E wind): Hs = {:.2} m, Tp = {:.2} s", 
             quebec_city.hs, quebec_city.tp);
    
    let tadoussac = st_lawrence_waves(25.0, 270.0, "Tadoussac");
    println!("Tadoussac (25 m/s W wind): Hs = {:.2} m, Tp = {:.2} s", 
             tadoussac.hs, tadoussac.tp);
    
    let sept_iles = st_lawrence_waves(22.0, 45.0, "Sept-Iles");
    println!("Sept-Iles (22 m/s NE wind): Hs = {:.2} m, Tp = {:.2} s", 
             sept_iles.hs, sept_iles.tp);
    
    // Demonstrate Wave Generation module
    println!("\n\n2. WAVE GENERATION MODULE");
    println!("=========================");
    
    // Pierson-Moskowitz calculations
    println!("\n2.1 Pierson-Moskowitz Spectrum:");
    let wind_speeds = [10.0, 15.0, 20.0, 25.0];
    for speed in wind_speeds {
        let hs = WaveGeneration::pierson_moskowitz_hs(speed);
        let fp = WaveGeneration::pierson_moskowitz_fp(speed);
        println!("Wind: {:.0} m/s → Hs = {:.2} m, fp = {:.4} Hz", speed, hs, fp);
    }
    
    // JONSWAP fetch-limited
    println!("\n2.2 JONSWAP Fetch-Limited Generation:");
    let conditions = WindConditions {
        speed: 15.0,
        direction: 270.0,
        fetch: 50_000.0,
        duration: 3600.0 * 4.0,
        water_depth: 100.0,
    };
    let jonswap_waves = WaveGeneration::jonswap_fetch_limited(&conditions);
    println!("Fetch-limited: Hs = {:.2} m, Tp = {:.2} s, Fully developed: {}", 
             jonswap_waves.hs, jonswap_waves.tp, jonswap_waves.is_fully_developed);
    
    // Duration-limited
    println!("\n2.3 Duration-Limited Generation:");
    let duration_waves = WaveGeneration::duration_limited(&conditions);
    println!("Duration-limited: Hs = {:.2} m, Tp = {:.2} s", 
             duration_waves.hs, duration_waves.tp);
    
    // Combined fetch-duration limited
    println!("\n2.4 Combined Fetch-Duration Limited:");
    let combined_waves = WaveGeneration::fetch_duration_limited(&conditions);
    println!("Combined: Hs = {:.2} m, Tp = {:.2} s", 
             combined_waves.hs, combined_waves.tp);
    
    // Wave dispersion calculations
    println!("\n2.5 Wave Dispersion:");
    let wavelengths = [50.0, 100.0, 200.0, 500.0];
    let depths = [10.0, 50.0, 200.0];
    
    for &wavelength in &wavelengths {
        let deep_speed = WaveDispersion::deep_Water_phase_speed(wavelength);
        println!("λ = {:.0} m:", wavelength);
        println!("  Deep water phase speed: {:.2} m/s", deep_speed);
        
        for &depth in &depths {
            let finite_speed = WaveDispersion::finite_depth_phase_speed(wavelength, depth);
            let group_speed = WaveDispersion::group_velocity(wavelength, depth);
            let is_deep = WaveDispersion::is_deep_water(wavelength, depth);
            let is_shallow = WaveDispersion::is_shallow_water(wavelength, depth);
            
            println!("  Depth {:.0} m: C = {:.2} m/s, Cg = {:.2} m/s (Deep: {}, Shallow: {})", 
                     depth, finite_speed, group_speed, is_deep, is_shallow);
        }
        println!();
    }
    
    // Wave packet propagation
    println!("2.6 Wave Packet Propagation:");
    let mut packet = WavePacket {
        position: (0.0, 0.0),
        energy: 1.0,
        frequency: 0.1,
        direction: 45.0,
        spreading: 15.0,
    };
    
    println!("Initial position: ({:.0}, {:.0})", packet.position.0, packet.position.1);
    for i in 1..=5 {
        packet.propagate(600.0, 50.0); // 10 minutes
        println!("After {} min: ({:.0}, {:.0}), Energy: {:.4}", 
                 i * 10, packet.position.0, packet.position.1, packet.energy);
    }
    
    // JONSWAP spectrum values
    println!("\n2.7 JONSWAP Spectrum:");
    let hs = 3.0;
    let tp = 8.0;
    let gamma = 3.3;
    let frequencies = [0.08, 0.10, 0.125, 0.15, 0.20]; // Hz
    
    for &f in &frequencies {
        let spectrum_value = jonswap_sepectrum(f, hs, tp, gamma);
        println!("f = {:.3} Hz: S(f) = {:.2} m²s", f, spectrum_value);
    }
    
    // Wave field demonstration
    println!("\n2.8 Wave Field Propagation:");
    let wind_field = WindConditions {
        speed: 20.0,
        direction: 270.0,
        fetch: 100_000.0,
        duration: 3600.0 * 6.0,
        water_depth: 50.0,
    };
    
    let bathymetry = Box::new(|x: f64, y: f64| {
        let r = (x * x + y * y).sqrt();
        50.0 + 30.0 * (-r / 50_000.0).exp()
    });
    
    let mut wave_field = WaveField::from_wind_field(
        &wind_field,
        (0.0, 0.0),
        (-50_000.0, 50_000.0, -50_000.0, 50_000.0),
        bathymetry,
    );
    
    println!("Initial wave field: {} wave packets", wave_field.packets.len());
    
    // Advance the field and sample statistics
    for hour in 1..=3 {
        wave_field.advance(3600.0); // 1 hour
        let stats = wave_field.get_wave_stats_at(25_000.0, 0.0, 5_000.0);
        match stats {
            Some((hs, tp)) => println!("Hour {}: Hs = {:.2} m, Tp = {:.2} s at (25km, 0)", hour, hs, tp),
            None => println!("Hour {}: No waves at sample point", hour),
        }
    }
    
    // Demonstrate visualization module
    println!("\n\n3. VISUALIZATION MODULE");
    println!("=======================");
    
    // Create visualization of the wave field
    match visualize_wave_field(&wave_field, "wave_field_demo.png") {
        Ok(()) => println!("Wave field visualization saved to wave_field_demo.png"),
        Err(e) => println!("Visualization error: {}", e),
    }
    
    // Create another wave field scenario for visualization
    let storm_conditions = WindConditions {
        speed: 30.0,
        direction: 45.0,
        fetch: 200_000.0,
        duration: 3600.0 * 12.0,
        water_depth: 100.0,
    };
    
    let storm_bathymetry = Box::new(|x: f64, y: f64| {
        let shore_dist = (x * x + y * y).sqrt();
        if shore_dist < 20_000.0 {
            10.0 + shore_dist / 2_000.0
        } else {
            100.0
        }
    });
    
    let mut storm_field = WaveField::from_wind_field(
        &storm_conditions,
        (-30_000.0, -30_000.0),
        (-100_000.0, 100_000.0, -100_000.0, 100_000.0),
        storm_bathymetry,
    );
    
    // Advance storm field
    storm_field.advance(3600.0 * 2.0); // 2 hours
    
    match visualize_wave_field(&storm_field, "storm_field_demo.png") {
        Ok(()) => println!("Storm wave field visualization saved to storm_field_demo.png"),
        Err(e) => println!("Storm visualization error: {}", e),
    }
    
    println!("\n=== All public functions demonstrated ===");
}
