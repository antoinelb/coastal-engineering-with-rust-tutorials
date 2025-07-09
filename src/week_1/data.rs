use crate::week_1::wave_stats::{WaveMeasurement, WaveTimeSeries};
use polars::io::ipc::{IpcReader, IpcWriter};
use polars::prelude::*;
use reqwest;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub async fn fetch_wave_data(
    token: &str,
    start_date: &str,
    end_date: &str,
) -> Result<(), Box<dyn Error>> {
    let raw_path = Path::new("data/ctd_data.json");
    let df_path = Path::new("data/ctd_data.ipc");

    download_ctd_data(raw_path, token, start_date, end_date).await?;
    read_ctd_data(raw_path, df_path)?;

    return Ok(());
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
    }

    return Ok(());
}

fn read_ctd_data(raw_path: &Path, df_path: &Path) -> Result<(), Box<dyn Error>> {
    if !df_path.exists() {
        let string_data = std::fs::read_to_string(raw_path).expect("Unable to read raw file.");
        let json_data: serde_json::Value =
            serde_json::from_str(&string_data).expect("Unable to parse raw data.");

        let mut dataframes: Vec<DataFrame> = Vec::new();

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

                        if datetimes.len() == vals.len() {
                            dataframes.push(df! {
                                "datetime" => datetimes,
                                sensor_name => vals,
                            }?)
                        }
                    }
                }
            }
        }

        if !dataframes.is_empty() {
            let mut dataframes_iter = dataframes.into_iter();
            let mut dataframe = dataframes_iter.next().unwrap();

            for dataframe_ in dataframes_iter.skip(1) {
                dataframe = dataframe
                    .join(
                        &dataframe_,
                        ["datetime"],
                        ["datetime"],
                        JoinArgs::new(JoinType::Inner),
                        None,
                    )
                    .unwrap();
            }

            let mut file = File::create(df_path)?;
            IpcWriter::new(&mut file).finish(&mut dataframe.clone())?;
        }
    }

    return Ok(());
}

fn calculate_water_level() {
    // Water_Level = (Pressure - Atmospheric_Pressure) / (Density × g)
}
