use plotters::prelude::*;

use crate::wave_generation::WaveField;

pub fn visualize_wave_field(
    field: &WaveField,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Wave Field Propagation", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(
            field.domain.0..field.domain.1,
            field.domain.2..field.domain.3,
        )?;

    chart
        .configure_mesh()
        .x_desc("Distance East (m)")
        .y_desc("Distance North (m)")
        .draw()?;

    // draw wave packets as circles with size proportional to energy
    for packet in &field.packets {
        let size = (packet.energy * 100.0).sqrt().max(1.0);
        chart.draw_series(std::iter::once(Circle::new(
            (packet.position.0, packet.position.1),
            size,
            ShapeStyle::from(BLUE).filled(),
        )))?;
    }

    root.present()?;

    return Ok(());
}
