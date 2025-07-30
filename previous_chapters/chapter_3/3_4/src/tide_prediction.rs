use chrono::{DateTime, Duration, Utc};
use std::f64::consts::PI;

// tidal constituent parameters
#[derive(Debug, Clone)]
pub struct TidalConstituent {
    pub name: String,
    pub frequency: f64, // cycles per hour
    pub amplitude: f64, // meters
    pub phase: f64,     // degrees
    pub speed: f64,     // degrees per hour
}

// collection of constituents for a location
#[derive(Debug, Clone)]
pub struct TidalStation {
    pub name: String,
    pub lat: f64,
    pub lon: f64,
    pub timezone: i32,
    pub datum_offset: f64, // mean sea level to chart datum
    pub constituents: Vec<TidalConstituent>,
}

#[derive(Debug, Clone)]
pub struct TideExtreme {
    pub time: DateTime<Utc>,
    pub height: f64,
    pub extreme_type: ExtremeType,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ExtremeType {
    High,
    Low,
}

impl TidalStation {
    pub fn predict_height(&self, time: DateTime<Utc>) -> Result<f64, String> {
        // reference time for phase calculation (usually start of year)
        let reference_time = DateTime::parse_from_rfc3339("2024-01-01T00:00:00z")
            .unwrap()
            .with_timezone(&Utc);

        let hours_since_reference = (time - reference_time).num_seconds() as f64 / 3600.0;

        // sum all constituents
        let mut height = self.datum_offset;

        for constituent in &self.constituents {
            let phase = constituent.phase * PI / 180.0; // convert to radians
            let omega = constituent.speed * PI / 180.0; // convert to radians per hour

            height += constituent.amplitude * (omega * hours_since_reference + phase).cos();
        }

        return Ok(height);
    }

    // find high and low tides in a time range
    pub fn find_extremes(
        &self,
        start_time: DateTime<Utc>,
        end_time: DateTime<Utc>,
    ) -> Result<Vec<TideExtreme>, String> {
        let mut extremes = Vec::new();
        let mut current_time = start_time;
        let time_step = Duration::minutes(6); // 6-minute intervals

        let previous_height = self.predict_height(current_time)?;
        current_time += time_step;
        let mut current_height = self.predict_height(current_time)?;
        let mut increasing = current_height > previous_height;

        while current_time < end_time {
            let next_time = current_time + time_step;
            let next_height = self.predict_height(next_time)?;

            let now_increasing = next_height > current_height;

            if increasing != now_increasing {
                // refine extremum time using quadratic interpolation
                let refined = self.refine_extremum(
                    current_time - time_step,
                    current_time,
                    current_time + time_step,
                )?;
                extremes.push(TideExtreme {
                    time: refined.0,
                    height: refined.1,
                    extreme_type: if increasing {
                        ExtremeType::High
                    } else {
                        ExtremeType::Low
                    },
                });
                increasing = now_increasing;
            }

            current_height = next_height;
            current_time = next_time;
        }

        return Ok(extremes);
    }

    fn refine_extremum(
        &self,
        t1: DateTime<Utc>,
        t2: DateTime<Utc>,
        t3: DateTime<Utc>,
    ) -> Result<(DateTime<Utc>, f64), String> {
        let h1 = self.predict_height(t1)?;
        let h2 = self.predict_height(t2)?;
        let h3 = self.predict_height(t3)?;

        // quadratic interpolation
        // fit parabola: h(t) = at^2 + bt + c
        let dt = (t2 - t1).num_seconds() as f64;
        let a = (h3 - 2.0 * h2 + h1) / (2.0 * dt * dt);
        let b = (h3 - h1) / (2.0 * dt);

        let t_offset = -b / (2.0 * a); // extremum at t = -b/(2a)
        let extremum_time = t2 + Duration::seconds(t_offset as i64);
        let extremum_height = self.predict_height(extremum_time)?;

        return Ok((extremum_time, extremum_height));
    }
}
