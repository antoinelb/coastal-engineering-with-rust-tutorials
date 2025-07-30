use crate::utils::print::{done_print, load_print};
use crate::wave_stats::{WaveMeasurement, WaveTimeSeries};
use chrono::{NaiveDate, NaiveDateTime};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;

#[derive(Deserialize, Serialize)]
struct CtdData {
    datetime: String,
    pressure: f64,
    sigma_t: f64,
}

pub async fn fetch_wave_data(
    token: &str,
    start_date: &str,
    end_date: &str,
) -> Result<WaveTimeSeries, Box<dyn Error>> {
    let raw_path = Path::new("data/ctd_data.json");
    let parsed_path = Path::new("data/parsed_ctd_data.json");
    let path = Path::new("data/data.json");

    let data: WaveTimeSeries;

    if !path.exists() {
        download_ctd_data(raw_path, token, start_date, end_date).await?;
        let ctd_data = parse_ctd_data(raw_path, parsed_path)?;
        // let mut ctd_data = read_ctd_data(raw_path, df_path)?;
        //
        data = calculate_water_level(ctd_data);

        let mut file = File::create(path)?;
        file.write_all(serde_json::to_string_pretty(&data).unwrap().as_bytes())?;
    } else {
        let string_data = std::fs::read_to_string(path)?;
        data = serde_json::from_str(&string_data)?;
    }

    return Ok(data);
}

async fn download_ctd_data(
    path: &Path,
    token: &str,
    start_date: &str,
    end_date: &str,
) -> Result<(), Box<dyn Error>> {
    let device_code = "CTD";
    let station_code = "SCVIP";

    if !path.exists() {
        load_print("Downloading CTD data...");

        let client = reqwest::Client::new();
        let url = "https://data.oceannetworks.ca/api/scalardata/location";

        let resp = client
            .get(url)
            .query(&[
                ("method", "get"),
                ("token", token),
                ("locationCode", station_code),
                ("deviceCategoryCode", device_code),
                ("dateFrom", start_date),
                ("dateTo", end_date),
            ])
            .send()
            .await?;

        let data: serde_json::Value = resp.json().await?;

        std::fs::create_dir_all(path.parent().unwrap()).unwrap();

        let mut file = File::create(path)?;
        file.write_all(serde_json::to_string_pretty(&data)?.as_bytes())?;

        done_print("Downloaded CTD data.");
    } else {
        done_print("CTD data already downloaded.");
    }

    return Ok(());
}

fn parse_ctd_data(raw_path: &Path, parsed_path: &Path) -> Result<Vec<CtdData>, Box<dyn Error>> {
    let mut data: Vec<CtdData>;

    if !parsed_path.exists() {
        load_print("Parsing CTD data...");
        let string_data = std::fs::read_to_string(raw_path).expect("Unable to read raw file.");
        let json_data: serde_json::Value =
            serde_json::from_str(&string_data).expect("Unable to parse raw data.");

        data = Vec::new();
        let mut pressures: HashMap<String, f64> = HashMap::new();
        let mut sigma_ts: HashMap<String, f64> = HashMap::new();

        if let Some(sensor_data) = json_data["sensorData"].as_array() {
            for sensor in sensor_data {
                if let (Some(sensor_name), Some(data)) =
                    (sensor["sensorName"].as_str(), sensor["data"].as_object())
                {
                    if let (Some(sample_times), Some(values)) =
                        (data["sampleTimes"].as_array(), data["values"].as_array())
                    {
                        let datetimes: Vec<String> = sample_times
                            .iter()
                            .filter_map(|t| t.as_str().map(|s| s.to_string()))
                            .collect();
                        let vals: Vec<f64> = values.iter().filter_map(|v| v.as_f64()).collect();

                        match sensor_name {
                            "Pressure" => {
                                pressures = HashMap::from_iter(
                                    datetimes
                                        .iter()
                                        .zip(vals.iter())
                                        .map(|(datetime, val)| (datetime.clone(), *val))
                                        .collect::<Vec<(String, f64)>>(),
                                );
                            }
                            "Sigma-t" => {
                                sigma_ts = HashMap::from_iter(
                                    datetimes
                                        .iter()
                                        .zip(vals.iter())
                                        .map(|(datetime, val)| (datetime.clone(), *val))
                                        .collect::<Vec<(String, f64)>>(),
                                )
                            }
                            _ => {}
                        }
                    }
                }
            }
        }

        for (datetime, pressure) in pressures.into_iter() {
            if let Some(&sigma_t) = sigma_ts.get(&datetime) {
                data.push(CtdData {
                    datetime,
                    pressure,
                    sigma_t,
                })
            }
        }

        let mut file = File::create(parsed_path)?;
        file.write_all(serde_json::to_string_pretty(&data).unwrap().as_bytes())?;

        done_print("Parsed CTD data.");
    } else {
        let string_data =
            std::fs::read_to_string(parsed_path).expect("Unable to read parsed file.");
        data = serde_json::from_str(&string_data)?;
        done_print("Read parsed CTD data.");
    }

    return Ok(data);
}

fn calculate_water_level(data: Vec<CtdData>) -> WaveTimeSeries {
    // Water_Level = (Pressure - Atmospheric_Pressure) / (Density × g)
    let measurements: Vec<WaveMeasurement> = data
        .iter()
        .map(|d| WaveMeasurement {
            time: (NaiveDateTime::parse_from_str(&d.datetime, "%Y-%m-%dT%H:%M:%S.%fZ").unwrap()
                - NaiveDate::from_ymd_opt(2024, 1, 1)
                    .unwrap()
                    .and_hms_opt(0, 0, 0)
                    .unwrap())
            .as_seconds_f64(),
            elevation: d.pressure - 101.325 / ((d.sigma_t + 1000.0) * 9.81),
        })
        .collect();
    let sampling_rate = measurements[1].time - measurements[0].time;

    return WaveTimeSeries {
        measurements,
        sampling_rate,
    };
}
