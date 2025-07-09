use crate::week_1::wave_stats::{Wave, WaveStatistics, WaveTimeSeries};
use plotters::prelude::*;

pub fn plot_wave_timeseries(
    data: &WaveTimeSeries,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let measurements: Vec<(f64, f64)> = data
        .measurements
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

    chart
        .configure_mesh()
        .x_desc("Time (s)")
        .y_desc("Elevation (m)")
        .draw()?;

    chart
        .draw_series(LineSeries::new(measurements, BLUE))?
        .label("Surface elevation")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], BLUE));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    return Ok(());
}

pub fn plot_wave_height_distribution(
    waves: &[Wave],
    stats: &WaveStatistics,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_height = waves.iter().map(|w| w.height).fold(0.0, f64::max);
    let bin_width = max_height / 20.0;
    let mut histogram = [0; 20];

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

    chart
        .configure_mesh()
        .x_desc("Wave Height (m)")
        .y_desc("Count")
        .draw()?;

    chart.draw_series(histogram.iter().enumerate().map(|(i, &count)| {
        Rectangle::new(
            [
                (i as f64 * bin_width, 0.0),
                ((i + 1) as f64 * bin_width, count as f64),
            ],
            BLUE.filled(),
        )
    }))?;

    chart.draw_series(vec![
        PathElement::new([(stats.hs, 0.0), (stats.hs, max_count)], RED),
        PathElement::new([(stats.hmax, 0.0), (stats.hmax, max_count)], GREEN),
    ])?;

    chart
        .draw_series(std::iter::once(PathElement::new([(0.0, 0.0)], RED)))?
        .label(format!("Hs = {:.2} m", stats.hs))
        .legend(|(x, y)| PathElement::new([(x, y), (x + 10, y)], RED));
    chart
        .draw_series(std::iter::once(PathElement::new([(0.0, 0.0)], GREEN)))?
        .label(format!("Hmax = {:.2} m", stats.hmax))
        .legend(|(x, y)| PathElement::new([(x, y), (x + 10, y)], GREEN));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    return Ok(());
}
