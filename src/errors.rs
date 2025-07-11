use std::error::Error;
use std::fmt;

#[derive(Debug)]
pub enum WaveError {
    InsufficientData,
    InvalidMeasurement(String),
}

impl fmt::Display for WaveError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            WaveError::InsufficientData => write!(f, "Insufficient wave data"),
            WaveError::InvalidMeasurement(msg) => write!(f, "Invalid measurement: {}", msg),
        }
    }
}

impl Error for WaveError {}
