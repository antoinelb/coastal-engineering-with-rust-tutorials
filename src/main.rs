mod consts;
mod data_generator;
mod errors;
mod spectral_analysis;
mod waves;

fn main() {
    println!("=== Coastal Dynamics Wave Analysis Demo ===\n");

    // Showcase data_generator functionality
    println!("1. Generating synthetic wave data...");
    let duration = 300.0; // 5 minutes
    let sampling_rate = 4.0; // 4 Hz
    let hs = 2.5; // significant wave height in meters
    let peak_period = 8.0; // peak period in seconds

    let wave_data = data_generator::generate_wave_data(duration, sampling_rate, hs, peak_period);
    println!("   Generated {} wave measurements over {} seconds", 
             wave_data.measurements.len(), duration);
    println!("   Sampling rate: {} Hz", sampling_rate);
    println!("   Target significant height: {:.2} m", hs);
    println!("   Target peak period: {:.2} s\n", peak_period);

    // Showcase spectral analysis functionality
    println!("2. Performing spectral analysis...");
    match wave_data.perform_spectral_analysis() {
        Ok(analysis) => {
            println!("   Peak frequency: {:.4} Hz", analysis.peak_frequency);
            println!("   Peak period: {:.2} s", analysis.peak_period);
            println!("   Frequency resolution: {:.4} Hz", 
                     analysis.frequencies[1] - analysis.frequencies[0]);
            println!("   Spectrum length: {} frequency bins", analysis.spectrum.len());

            // Demonstrate sea-swell separation
            println!("\n3. Separating sea and swell components...");
            let wind_speed = 15.0; // m/s
            let components = analysis.separate_sea_swell(wind_speed);
            println!("   Wind speed: {:.1} m/s", wind_speed);
            println!("   Separation frequency: {:.4} Hz", components.separation_frequency);
            println!("   Sea Hs: {:.2} m", components.sea_hs);
            println!("   Swell Hs: {:.2} m", components.swell_hs);

            // Calculate spectral moments
            println!("\n4. Computing spectral moments...");
            if let Ok(m0) = wave_data.spectral_moments(0) {
                println!("   0th moment (variance): {:.4} m²", m0);
                println!("   Hs from spectrum: {:.2} m", 4.0 * m0.sqrt());
            }
            if let Ok(m1) = wave_data.spectral_moments(1) {
                println!("   1st moment: {:.4} Hz·m²", m1);
            }
            if let Ok(m2) = wave_data.spectral_moments(2) {
                println!("   2nd moment: {:.4} Hz²·m²", m2);
            }
        }
        Err(e) => {
            println!("   Error in spectral analysis: {:?}", e);
        }
    }

    println!("\n=== Demo Complete ===");
}
