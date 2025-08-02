use crate::wave_climate::WaveObservation;
use chrono::{DateTime, Utc};
use serde_json::Value;

pub struct OncDataFetcher {
    api_token: String,
    base_url: String,
}

impl OncDataFetcher {
    pub fn new(api_token: String) -> Self {
        OncDataFetcher {
            api_token,
            base_url: "https://data.oceannetworks.ca/api".to_string(),
        }
    }

    // fetch wave data for a specific location and time range
    pub async fn fetch_wave_data(
        &self,
        location_code: &str,
        start: DateTime<Utc>,
        end: DateTime<Utc>,
    ) -> Result<Vec<WaveObservation>, Box<dyn std::error::Error>> {
        let url = format!(
            "{}/scalardata/location/{}/devicecategory/WAVEMONITOR",
            self.base_url, location_code
        );

        let params = [
            ("token", &self.api_token),
            (
                "dateFrom",
                &start.format("%Y-%m-%dT%H:%M:%S.000Z").to_string(),
            ),
            ("dateTo", &end.format("%Y-%m-%dT%H:%M:%S.000Z").to_string()),
        ];

        let resp = reqwest::Client::new()
            .get(&url)
            .query(&params)
            .send()
            .await?
            .json::<Value>()
            .await?;

        let mut observations = Vec::new();

        if let Some(data) = resp["data"].as_array() {
            for record in data {
                if let (Some(timestamp), Some(hs), Some(tp), Some(dp)) = (
                    record["time"].as_str(),
                    record["significantWaveHeight"].as_f64(),
                    record["peakPeriod"].as_f64(),
                    record["peakDirection"].as_f64(),
                ) {
                    let timestamp = DateTime::parse_from_rfc3339(timestamp)?.with_timezone(&Utc);
                    observations.push(WaveObservation {
                        timestamp,
                        hs,
                        tp,
                        tm: record["meanPeriod"].as_f64().unwrap_or(tp * 0.8),
                        dp,
                        dm: record["meanDirection"].as_f64().unwrap_or(dp),
                        spread: record["directionalSpread"].as_f64().unwrap_or(30.0),
                        wind_speed: record["windSpeed"].as_f64(),
                        wind_dir: record["windDirection"].as_f64(),
                    });
                }
            }
        }

        Ok(observations)
    }
}
