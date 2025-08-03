use ndarray::{s, Array1, Array2};
use ndarray_linalg::{Inverse, Solve};
use serde::{Deserialize, Serialize};
use statrs::distribution::{Continuous, ContinuousCDF};
use std::collections::HashMap;
use std::f64::consts::{E, PI};

// GEV distribution parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GevParams {
    pub location: f64, // μ
    pub scale: f64,    // σ
    pub shape: f64,    // ξ
}

// GPD parameters for POT analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GpdParams {
    pub threshold: f64, // μ
    pub scale: f64,     // σ
    pub shape: f64,     // ξ
    pub rate: f64,      // λ (exceedances per year)
}

// confidence intervals for parameters and return values
#[derive(Debug, Clone)]
pub struct ConfidenceInterval {
    pub estimate: f64,
    pub lower: f64, // e.g. 5th percentile
    pub upper: f64, // e.g. 95th percentile
    pub std_error: f64,
}

#[derive(Debug, Clone)]
pub enum ExtremeValueParams {
    Gev(GevParams),
    Gpd(GpdParams),
}

// results from extreme value analysis
#[derive(Debug)]
pub struct EvaResults {
    pub params: ExtremeValueParams,
    pub return_values: HashMap<u32, ConfidenceInterval>, // return period -> estimate
}
