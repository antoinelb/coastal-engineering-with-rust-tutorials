use crate::wave_climate::{DirectionalBin, WaveClimateAnalyzer};
use plotters::{prelude::*, style::full_palette::ORANGE};

pub fn plot_wave_rose(
    bins: &[DirectionalBin],
    output_path: &str,
    title: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (800, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_radius = bins
        .iter()
        .map(|b| b.hs_mean * b.occurrence / 100.0)
        .fold(0.0, f64::max)
        * 1.2;

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 40))
        .margin(20)
        .build_cartesian_2d(-max_radius..max_radius, -max_radius..max_radius)?;

    for bin in bins {
        let theta = (90.0 - bin.direction).to_radians();
        let r = bin.hs_mean * bin.occurrence / 100.0;

        // create wedge polygon
        let n_points = 20;
        let dtheta = bin.width.to_radians() / n_points as f64;
        let start_theta = theta - bin.width.to_radians() / 2.0;

        let mut points = vec![(0.0, 0.0)];
        for i in 0..=n_points {
            let angle = start_theta + i as f64 * dtheta;
            points.push((r * angle.cos(), r * angle.sin()));
        }

        // color based on mean hs
        let color = match bin.hs_mean {
            h if h < 1.0 => BLUE.mix(0.3),
            h if h < 2.0 => BLUE.mix(0.5),
            h if h < 3.0 => YELLOW.mix(0.7),
            h if h < 4.0 => ORANGE.mix(0.7),
            _ => RED.mix(0.7),
        };

        chart.draw_series(std::iter::once(Polygon::new(points, color)))?;
    }

    // add compas directions
    let directions: Vec<(&str, f64)> = vec![("N", 0.0), ("E", 90.0), ("S", 180.0), ("W", 270.0)];
    for (label, deg) in directions {
        let theta = (90.0 - deg).to_radians();
        let x = max_radius * 0.9 * theta.cos();
        let y = max_radius * 0.9 * theta.sin();

        chart.draw_series(std::iter::once(Text::new(
            label,
            (x, y),
            ("sans-serif", 20).into_font(),
        )))?;
    }

    root.present()?;
    Ok(())
}

pub fn plot_seasonal_comparison(
    analyzer: &WaveClimateAnalyzer,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let seasonal = analyzer.seasonal_analysis();
    let root = BitMapBackend::new(output_path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Seasonal Wave Climate - BC Coast", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(["Winter", "Spring", "Summer", "Autumn"].as_ref(), 0.0..6.0)?;

    chart
        .configure_mesh()
        .y_desc("Significant Wave Height (m)")
        .x_desc("Season")
        .draw()?;

    let seasons = [
        ("Winter", &seasonal.winter),
        ("Spring", &seasonal.spring),
        ("Summer", &seasonal.summer),
        ("Autumn", &seasonal.autumn),
    ];

    chart
        .draw_series(seasons.iter().map(|(name, stats)| {
            Rectangle::new([(name, 0.0), (name, stats.mean_hs)], BLUE.filled())
        }))?
        .label("Mean Hs")
        .legend(|(x, y)| Rectangle::new([(x, y), (x + 10, y + 10)], BLUE.filled()));

    chart
        .draw_series(seasons.iter().map(|(name, stats)| {
            let p90 = stats.percentiles.get(&90).unwrap_or(&0.0);
            Circle::new((name, *p90), 5, RED.filled())
        }))?
        .label("90th percentile")
        .legend(|(x, y)| Circle::new((x + 5, y + 5), 5, RED.filled()));

    chart
        .draw_series(seasons.iter().map(|(name, stats)| {
            PathElement::new(vec![(name, stats.max_hs)], GREEN.stroke_width(2))
        }))?
        .label("Maximum")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], GREEN));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;
    Ok(())
}
