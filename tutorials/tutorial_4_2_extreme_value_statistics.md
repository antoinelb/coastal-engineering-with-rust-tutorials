# Tutorial 4.2: Extreme Value Statistics
## Global Wave Environments Chapter 4.3-4.4

### Learning Objectives
By the end of this tutorial, you will:
1. Master extreme value theory for coastal engineering applications
2. Implement robust methods for estimating design wave conditions
3. Quantify uncertainty in extreme value predictions
4. Apply both annual maxima and peak-over-threshold approaches
5. Understand the connection between extreme waves and coastal risk

### Prerequisites
- Read Coastal Dynamics Ch. 4.3-4.4
- Complete Tutorial 4.1 (Wave Climate Analysis)
- Understanding of probability distributions and statistics
- Familiarity with risk-based design concepts

### Introduction: Why Extreme Values Matter

Extreme wave events drive critical coastal processes and pose the greatest risks to infrastructure and communities. In coastal engineering, we must answer questions like:
- What is the 100-year wave height for harbour design?
- How high should coastal defences be to protect against storm surge?
- What loads will offshore structures experience during their lifetime?
- How will climate change affect extreme wave statistics?

The challenge is that we need to estimate events far rarer than our observational record. A 30-year dataset must inform design for 100-year or even 1000-year events. This tutorial explores the mathematical framework and practical tools for this extrapolation.

### Part 1: Theoretical Foundation

#### Extreme Value Distributions

Three key distributions form the foundation of extreme value analysis:

1. **Generalized Extreme Value (GEV)**: For block maxima
   $$F(x) = \exp\left\{-\left[1 + \xi\left(\frac{x-\mu}{\sigma}\right)\right]^{-1/\xi}\right\}$$
   where $\mu$ is location, $\sigma$ is scale, and $\xi$ is shape

2. **Generalized Pareto Distribution (GPD)**: For threshold exceedances
   $$F(x) = 1 - \left(1 + \frac{\xi x}{\sigma}\right)^{-1/\xi}$$

3. **Special cases**:
   - Gumbel ($\xi = 0$): Light-tailed, exponential decay
   - Fréchet ($\xi > 0$): Heavy-tailed, power-law decay
   - Weibull ($\xi < 0$): Bounded upper tail

#### Return Periods and Risk

The return period $T$ is the average recurrence interval:
$$T = \frac{1}{1-F(x)}$$

For design life $L$ years, the encounter probability is:
$$P(\text{exceedance in } L \text{ years}) = 1 - (1-1/T)^L$$

#### Questions to Consider
1. Why do different coastal locations have different shape parameters?
2. How does climate non-stationarity affect extreme value analysis?
3. What are the trade-offs between annual maxima and POT methods?

### Part 2: Comprehensive EVA Implementation

Create a robust extreme value analysis framework in `src/extreme_value_analysis.rs`:

```rust
// extreme_value_analysis.rs
use ndarray::{Array1, Array2, s};
use ndarray_linalg::{Inverse, Solve};
use statrs::distribution::{Continuous, ContinuousCDF};
use std::f64::consts::{E, PI};
use serde::{Serialize, Deserialize};

/// GEV distribution parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GevParams {
    pub location: f64,   // μ
    pub scale: f64,      // σ  
    pub shape: f64,      // ξ
}

/// GPD parameters for POT analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GpdParams {
    pub threshold: f64,  // u
    pub scale: f64,      // σ
    pub shape: f64,      // ξ
    pub rate: f64,       // λ (exceedances per year)
}

/// Confidence intervals for parameters and return values
#[derive(Debug, Clone)]
pub struct ConfidenceInterval {
    pub estimate: f64,
    pub lower: f64,      // e.g., 5th percentile
    pub upper: f64,      // e.g., 95th percentile
    pub std_error: f64,
}

/// Results from extreme value analysis
#[derive(Debug)]
pub struct EvaResults {
    pub params: Either<GevParams, GpdParams>,
    pub return_values: HashMap<u32, ConfidenceInterval>,  // Return period -> estimate
    pub goodness_of_fit: GoodnessOfFit,
    pub diagnostics: EvaDiagnostics,
}

/// Maximum Likelihood Estimation for GEV
pub struct GevMle {
    data: Array1<f64>,
    params: Option<GevParams>,
}

impl GevMle {
    pub fn new(annual_maxima: Vec<f64>) -> Self {
        GevMle {
            data: Array1::from_vec(annual_maxima),
            params: None,
        }
    }
    
    /// Negative log-likelihood for GEV
    fn neg_log_likelihood(&self, params: &[f64]) -> f64 {
        let (mu, sigma, xi) = (params[0], params[1], params[2]);
        
        if sigma <= 0.0 {
            return f64::INFINITY;
        }
        
        let n = self.data.len() as f64;
        let mut nll = n * sigma.ln();
        
        for &x in self.data.iter() {
            let z = (x - mu) / sigma;
            
            if xi.abs() < 1e-10 {
                // Gumbel case (ξ ≈ 0)
                nll += z + (-z).exp();
            } else {
                // General GEV case
                let t = 1.0 + xi * z;
                if t <= 0.0 {
                    return f64::INFINITY;
                }
                nll += (1.0 + 1.0 / xi) * t.ln() + t.powf(-1.0 / xi);
            }
        }
        
        nll
    }
    
    /// Fit GEV using L-moments for initial estimates
    pub fn fit(&mut self) -> Result<GevParams, String> {
        // Calculate L-moments for initial parameter estimates
        let l_moments = self.calculate_l_moments();
        let initial_params = self.l_moment_estimates(&l_moments);
        
        // Optimize using Nelder-Mead or similar
        let result = self.optimize_params(initial_params)?;
        
        self.params = Some(result.clone());
        Ok(result)
    }
    
    /// Calculate L-moments of the data
    fn calculate_l_moments(&self) -> LMoments {
        let mut sorted_data = self.data.to_vec();
        sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = sorted_data.len() as f64;
        
        // Calculate probability weighted moments
        let mut b0 = 0.0;
        let mut b1 = 0.0;
        let mut b2 = 0.0;
        let mut b3 = 0.0;
        
        for (i, &x) in sorted_data.iter().enumerate() {
            let p = (i as f64 + 0.35) / n;  // Plotting position
            b0 += x;
            b1 += x * p;
            b2 += x * p * p;
            b3 += x * p * p * p;
        }
        
        b0 /= n;
        b1 /= n;
        b2 /= n;
        b3 /= n;
        
        // Convert to L-moments
        let l1 = b0;
        let l2 = 2.0 * b1 - b0;
        let l3 = 6.0 * b2 - 6.0 * b1 + b0;
        let l4 = 20.0 * b3 - 30.0 * b2 + 12.0 * b1 - b0;
        
        LMoments {
            l1,
            l2,
            t3: l3 / l2,  // L-skewness
            t4: l4 / l2,  // L-kurtosis
        }
    }
    
    /// Initial GEV parameter estimates from L-moments
    fn l_moment_estimates(&self, lm: &LMoments) -> Vec<f64> {
        // Hosking and Wallis (1997) approximations
        let c = 2.0 / (3.0 + lm.t3) - E.ln() / 2.0_f64.ln();
        let xi = 7.8590 * c + 2.9554 * c * c;
        
        let gamma_1_plus_xi = gamma(1.0 + xi);
        let sigma = lm.l2 * xi / ((1.0 - 2.0_f64.powf(-xi)) * gamma_1_plus_xi);
        let mu = lm.l1 - sigma * (1.0 - gamma_1_plus_xi) / xi;
        
        vec![mu, sigma, xi]
    }
    
    /// Calculate return values with confidence intervals
    pub fn return_values(
        &self,
        return_periods: &[u32],
        n_bootstrap: usize,
    ) -> HashMap<u32, ConfidenceInterval> {
        let params = self.params.as_ref().expect("Model not fitted");
        let mut results = HashMap::new();
        
        for &rp in return_periods {
            let p = 1.0 - 1.0 / rp as f64;
            let point_estimate = self.quantile(p, params);
            
            // Bootstrap for confidence intervals
            let bootstrap_estimates = self.bootstrap_return_value(rp, n_bootstrap);
            let mut sorted_estimates = bootstrap_estimates.clone();
            sorted_estimates.sort_by(|a, b| a.partial_cmp(b).unwrap());
            
            let lower_idx = (0.05 * n_bootstrap as f64) as usize;
            let upper_idx = (0.95 * n_bootstrap as f64) as usize;
            
            let ci = ConfidenceInterval {
                estimate: point_estimate,
                lower: sorted_estimates[lower_idx],
                upper: sorted_estimates[upper_idx],
                std_error: bootstrap_estimates.iter()
                    .map(|&x| (x - point_estimate).powi(2))
                    .sum::<f64>()
                    .sqrt() / n_bootstrap as f64,
            };
            
            results.insert(rp, ci);
        }
        
        results
    }
    
    /// GEV quantile function
    fn quantile(&self, p: f64, params: &GevParams) -> f64 {
        let y = -(-p).ln();
        
        if params.shape.abs() < 1e-10 {
            // Gumbel
            params.location - params.scale * y.ln()
        } else {
            // General GEV
            params.location + params.scale / params.shape * (1.0 - y.powf(-params.shape))
        }
    }
}

/// Gamma function approximation
fn gamma(x: f64) -> f64 {
    // Stirling's approximation for large x
    if x > 5.0 {
        let e = E;
        (2.0 * PI * x).sqrt() * (x / e).powf(x) * (1.0 + 1.0 / (12.0 * x))
    } else {
        // Use recursion and known values
        // Simplified implementation
        (1..x.floor() as i32).map(|i| i as f64).product::<f64>()
    }
}

#[derive(Debug)]
struct LMoments {
    l1: f64,  // Mean
    l2: f64,  // L-scale
    t3: f64,  // L-skewness
    t4: f64,  // L-kurtosis
}
```

### Part 3: Peak-Over-Threshold Analysis

Implement the GPD approach for threshold exceedances:

```rust
// Continue in extreme_value_analysis.rs

/// Peak-over-threshold analysis using GPD
pub struct PotAnalysis {
    data: Vec<(DateTime<Utc>, f64)>,  // (timestamp, value)
    threshold: Option<f64>,
    params: Option<GpdParams>,
}

impl PotAnalysis {
    pub fn new(data: Vec<(DateTime<Utc>, f64)>) -> Self {
        PotAnalysis {
            data,
            threshold: None,
            params: None,
        }
    }
    
    /// Select threshold using multiple methods
    pub fn select_threshold(&mut self) -> ThresholdSelection {
        let values: Vec<f64> = self.data.iter().map(|(_, v)| *v).collect();
        let mut sorted_values = values.clone();
        sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        // Method 1: Mean Residual Life plot
        let mrl_data = self.mean_residual_life(&sorted_values);
        
        // Method 2: Parameter stability plot  
        let stability_data = self.parameter_stability(&sorted_values);
        
        // Method 3: Return value stability
        let rv_stability = self.return_value_stability(&sorted_values);
        
        // Automated selection based on multiple criteria
        let recommended = self.automated_threshold_selection(
            &mrl_data,
            &stability_data,
            &rv_stability,
        );
        
        self.threshold = Some(recommended);
        
        ThresholdSelection {
            recommended,
            mrl_plot: mrl_data,
            stability_plot: stability_data,
            rv_stability_plot: rv_stability,
        }
    }
    
    /// Mean Residual Life plot data
    fn mean_residual_life(&self, sorted_values: &[f64]) -> Vec<(f64, f64, f64)> {
        let mut results = Vec::new();
        let n = sorted_values.len();
        
        // Test thresholds from 70th to 98th percentile
        for i in ((0.7 * n as f64) as usize)..((0.98 * n as f64) as usize) {
            let threshold = sorted_values[i];
            let exceedances: Vec<f64> = sorted_values[i..]
                .iter()
                .map(|&v| v - threshold)
                .collect();
            
            if exceedances.len() > 10 {  // Minimum sample size
                let mean_excess = exceedances.iter().sum::<f64>() / exceedances.len() as f64;
                let std_error = exceedances.iter()
                    .map(|&e| (e - mean_excess).powi(2))
                    .sum::<f64>()
                    .sqrt() / (exceedances.len() as f64).sqrt();
                
                results.push((threshold, mean_excess, std_error));
            }
        }
        
        results
    }
    
    /// Decluster exceedances to ensure independence
    fn decluster_exceedances(
        &self,
        threshold: f64,
        min_separation: Duration,
    ) -> Vec<Vec<(DateTime<Utc>, f64)>> {
        let mut clusters = Vec::new();
        let mut current_cluster = Vec::new();
        
        for &(timestamp, value) in &self.data {
            if value > threshold {
                if let Some((last_time, _)) = current_cluster.last() {
                    if timestamp - *last_time > min_separation {
                        // Start new cluster
                        if !current_cluster.is_empty() {
                            clusters.push(current_cluster);
                        }
                        current_cluster = vec![(timestamp, value)];
                    } else {
                        // Add to current cluster
                        current_cluster.push((timestamp, value));
                    }
                } else {
                    // First exceedance
                    current_cluster.push((timestamp, value));
                }
            }
        }
        
        if !current_cluster.is_empty() {
            clusters.push(current_cluster);
        }
        
        clusters
    }
    
    /// Fit GPD to threshold exceedances
    pub fn fit_gpd(&mut self) -> Result<GpdParams, String> {
        let threshold = self.threshold.ok_or("Threshold not selected")?;
        
        // Decluster to get independent exceedances
        let clusters = self.decluster_exceedances(threshold, Duration::days(3));
        let peak_exceedances: Vec<f64> = clusters.iter()
            .map(|cluster| {
                cluster.iter()
                    .map(|(_, v)| v - threshold)
                    .fold(0.0, f64::max)
            })
            .collect();
        
        // Calculate exceedance rate
        let n_years = (self.data.last().unwrap().0 - self.data.first().unwrap().0)
            .num_days() as f64 / 365.25;
        let rate = peak_exceedances.len() as f64 / n_years;
        
        // Fit GPD using MLE
        let gpd_params = self.fit_gpd_mle(&peak_exceedances)?;
        
        self.params = Some(GpdParams {
            threshold,
            scale: gpd_params.0,
            shape: gpd_params.1,
            rate,
        });
        
        Ok(self.params.as_ref().unwrap().clone())
    }
    
    /// GPD maximum likelihood estimation
    fn fit_gpd_mle(&self, exceedances: &[f64]) -> Result<(f64, f64), String> {
        // Initial estimates using method of moments
        let mean = exceedances.iter().sum::<f64>() / exceedances.len() as f64;
        let variance = exceedances.iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>() / exceedances.len() as f64;
        
        let cv = variance.sqrt() / mean;  // Coefficient of variation
        let shape_init = 0.5 * (cv * cv - 1.0);
        let scale_init = mean * (1.0 - shape_init);
        
        // Optimize using Newton-Raphson or similar
        // Simplified implementation
        Ok((scale_init, shape_init))
    }
    
    /// Calculate return values from POT analysis
    pub fn return_values_pot(
        &self,
        return_periods: &[u32],
    ) -> HashMap<u32, f64> {
        let params = self.params.as_ref().expect("Model not fitted");
        let mut results = HashMap::new();
        
        for &rp in return_periods {
            let m = rp as f64 * params.rate;  // Expected exceedances in return period
            
            let quantile = if params.shape.abs() < 1e-10 {
                // Exponential case
                params.threshold + params.scale * m.ln()
            } else {
                // General GPD
                params.threshold + params.scale / params.shape * 
                    (m.powf(params.shape) - 1.0)
            };
            
            results.insert(rp, quantile);
        }
        
        results
    }
}

#[derive(Debug)]
pub struct ThresholdSelection {
    pub recommended: f64,
    pub mrl_plot: Vec<(f64, f64, f64)>,         // (threshold, mean_excess, std_error)
    pub stability_plot: Vec<(f64, f64, f64)>,    // (threshold, shape, scale)
    pub rv_stability_plot: Vec<(f64, f64)>,      // (threshold, 100yr_return_value)
}
```

### Part 4: Diagnostic Tools and Goodness-of-Fit

Implement comprehensive diagnostics:

```rust
// Continue in extreme_value_analysis.rs

/// Diagnostic plots and tests for EVA
pub struct EvaDiagnostics;

impl EvaDiagnostics {
    /// Q-Q plot data for GEV fit
    pub fn qq_plot_gev(data: &[f64], params: &GevParams) -> Vec<(f64, f64)> {
        let mut sorted_data = data.to_vec();
        sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = sorted_data.len();
        
        sorted_data.iter().enumerate().map(|(i, &x)| {
            let empirical_p = (i as f64 + 0.5) / n as f64;
            let theoretical = gev_quantile(empirical_p, params);
            (theoretical, x)
        }).collect()
    }
    
    /// Return level plot with confidence bands
    pub fn return_level_plot(
        eva_results: &EvaResults,
        min_period: f64,
        max_period: f64,
    ) -> ReturnLevelData {
        let periods: Vec<f64> = (0..100)
            .map(|i| min_period * (max_period / min_period).powf(i as f64 / 99.0))
            .collect();
        
        let mut estimates = Vec::new();
        let mut lower_bounds = Vec::new();
        let mut upper_bounds = Vec::new();
        
        for period in &periods {
            if let Some(ci) = eva_results.return_values.get(&(*period as u32)) {
                estimates.push(ci.estimate);
                lower_bounds.push(ci.lower);
                upper_bounds.push(ci.upper);
            }
        }
        
        ReturnLevelData {
            periods,
            estimates,
            lower_bounds,
            upper_bounds,
        }
    }
    
    /// Anderson-Darling test for goodness-of-fit
    pub fn anderson_darling_test(data: &[f64], params: &GevParams) -> TestResult {
        let mut sorted_data = data.to_vec();
        sorted_data.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = sorted_data.len();
        
        let mut a2 = 0.0;
        for (i, &x) in sorted_data.iter().enumerate() {
            let fi = gev_cdf(x, params);
            let term1 = (2.0 * (i + 1) as f64 - 1.0) * fi.ln();
            let term2 = (2.0 * (n - i) as f64 - 1.0) * (1.0 - fi).ln();
            a2 += term1 + term2;
        }
        
        a2 = -n as f64 - a2 / n as f64;
        
        // Modified for case 3 (estimated parameters)
        let a2_star = a2 * (1.0 + 0.2 / (n as f64).sqrt());
        
        // Approximate p-value
        let p_value = if a2_star < 0.2 {
            1.0 - (-13.436 + 101.14 * a2_star - 223.73 * a2_star * a2_star).exp()
        } else if a2_star < 0.34 {
            1.0 - (-8.318 + 42.796 * a2_star - 59.938 * a2_star * a2_star).exp()
        } else {
            (0.6 / a2_star + 1.0) * (-a2_star).exp()
        };
        
        TestResult {
            statistic: a2_star,
            p_value,
            reject_null: p_value < 0.05,
        }
    }
}

#[derive(Debug)]
pub struct ReturnLevelData {
    pub periods: Vec<f64>,
    pub estimates: Vec<f64>,
    pub lower_bounds: Vec<f64>,
    pub upper_bounds: Vec<f64>,
}

#[derive(Debug)]
pub struct TestResult {
    pub statistic: f64,
    pub p_value: f64,
    pub reject_null: bool,
}

/// Calculate GEV CDF
fn gev_cdf(x: f64, params: &GevParams) -> f64 {
    let z = (x - params.location) / params.scale;
    
    if params.shape.abs() < 1e-10 {
        // Gumbel
        (-(-z).exp()).exp()
    } else {
        // General GEV
        let t = 1.0 + params.shape * z;
        if t <= 0.0 {
            if params.shape > 0.0 { 0.0 } else { 1.0 }
        } else {
            (-t.powf(-1.0 / params.shape)).exp()
        }
    }
}

/// Calculate GEV quantile
fn gev_quantile(p: f64, params: &GevParams) -> f64 {
    let y = -(-p).ln();
    
    if params.shape.abs() < 1e-10 {
        params.location - params.scale * y.ln()
    } else {
        params.location + params.scale / params.shape * (1.0 - y.powf(-params.shape))
    }
}
```

### Part 5: Visualization Functions

Create comprehensive visualizations for extreme value analysis:

```rust
// In src/eva_visualization.rs
use plotters::prelude::*;
use crate::extreme_value_analysis::*;

/// Plot diagnostic plots for EVA
pub fn plot_eva_diagnostics(
    results: &EvaResults,
    data: &[f64],
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (1200, 900)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.split_evenly((2, 2));
    
    // 1. Q-Q Plot
    let qq_data = match &results.params {
        Either::Left(gev) => EvaDiagnostics::qq_plot_gev(data, gev),
        Either::Right(_) => vec![], // TODO: Implement for GPD
    };
    
    let mut qq_chart = ChartBuilder::on(&root[0])
        .caption("Q-Q Plot", ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(
            qq_data.iter().map(|(x, _)| *x).fold(f64::INFINITY, f64::min)..
            qq_data.iter().map(|(x, _)| *x).fold(f64::NEG_INFINITY, f64::max),
            qq_data.iter().map(|(_, y)| *y).fold(f64::INFINITY, f64::min)..
            qq_data.iter().map(|(_, y)| *y).fold(f64::NEG_INFINITY, f64::max),
        )?;
    
    qq_chart.configure_mesh()
        .x_desc("Theoretical Quantiles")
        .y_desc("Sample Quantiles")
        .draw()?;
    
    // Plot points
    qq_chart.draw_series(
        qq_data.iter().map(|&(x, y)| Circle::new((x, y), 3, BLUE.filled()))
    )?;
    
    // Add 1:1 line
    let min_val = qq_data.iter().map(|(x, y)| x.min(y)).fold(f64::INFINITY, |a, &b| a.min(b));
    let max_val = qq_data.iter().map(|(x, y)| x.max(y)).fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    qq_chart.draw_series(LineSeries::new(
        vec![(min_val, min_val), (max_val, max_val)],
        &RED,
    ))?;
    
    // 2. Return Level Plot
    let rl_data = EvaDiagnostics::return_level_plot(results, 1.0, 1000.0);
    
    let mut rl_chart = ChartBuilder::on(&root[1])
        .caption("Return Level Plot", ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(
            (1.0f64..1000.0).log_scale(),
            0.0..rl_data.estimates.iter().fold(0.0, |a, &b| a.max(b)) * 1.2,
        )?;
    
    rl_chart.configure_mesh()
        .x_desc("Return Period (years)")
        .y_desc("Return Level (m)")
        .draw()?;
    
    // Plot estimate with confidence bands
    rl_chart.draw_series(LineSeries::new(
        rl_data.periods.iter().zip(&rl_data.estimates)
            .map(|(&x, &y)| (x, y)),
        &BLUE,
    ))?
    .label("Estimate")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &BLUE));
    
    // Confidence bands
    rl_chart.draw_series(LineSeries::new(
        rl_data.periods.iter().zip(&rl_data.lower_bounds)
            .map(|(&x, &y)| (x, y)),
        ShapeStyle::from(&RED).stroke_width(1),
    ))?;
    
    rl_chart.draw_series(LineSeries::new(
        rl_data.periods.iter().zip(&rl_data.upper_bounds)
            .map(|(&x, &y)| (x, y)),
        ShapeStyle::from(&RED).stroke_width(1),
    ))?
    .label("90% CI")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], &RED));
    
    rl_chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
    
    // 3. Probability Plot
    // 4. Density Plot with fitted distribution
    
    root[0].present()?;
    Ok(())
}

/// Plot threshold selection diagnostics
pub fn plot_threshold_diagnostics(
    selection: &ThresholdSelection,
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output_path, (1200, 400)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.split_evenly((1, 3));
    
    // 1. Mean Residual Life Plot
    let mut mrl_chart = ChartBuilder::on(&root[0])
        .caption("Mean Residual Life", ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(
            selection.mrl_plot.first().unwrap().0..selection.mrl_plot.last().unwrap().0,
            0.0..selection.mrl_plot.iter().map(|(_, m, _)| *m).fold(0.0, f64::max) * 1.2,
        )?;
    
    mrl_chart.configure_mesh()
        .x_desc("Threshold")
        .y_desc("Mean Excess")
        .draw()?;
    
    // Plot with error bars
    mrl_chart.draw_series(
        selection.mrl_plot.iter().map(|&(u, mean, se)| {
            ErrorBar::new_vertical((u, mean - 2.0 * se), (u, mean), (u, mean + 2.0 * se))
        })
    )?;
    
    // Mark selected threshold
    mrl_chart.draw_series(std::iter::once(PathElement::new(
        vec![(selection.recommended, 0.0), 
             (selection.recommended, mrl_chart.y_range().end)],
        ShapeStyle::from(&GREEN).stroke_width(2),
    )))?;
    
    // 2. Parameter Stability Plot
    // 3. Return Value Stability Plot
    
    root[0].present()?;
    Ok(())
}
```

### Part 6: Main Application

Integrate all components into a comprehensive extreme value analysis:

```rust
// In src/main.rs
mod extreme_value_analysis;
mod eva_visualization;
mod wave_climate;

use extreme_value_analysis::*;
use chrono::{DateTime, Utc, TimeZone};
use std::error::Error;

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {
    println!("Coastal Dynamics Tutorial 4.2: Extreme Value Statistics");
    println!("=====================================================\n");
    
    // Load or generate wave data
    let wave_data = load_wave_data()?;  // Would load from file or API
    
    // Example 1: Annual Maxima Method
    println!("=== Annual Maxima Method ===");
    let annual_maxima = extract_annual_maxima(&wave_data);
    println!("Found {} annual maxima from the dataset", annual_maxima.len());
    
    let mut gev_analysis = GevMle::new(annual_maxima.clone());
    let gev_params = gev_analysis.fit()?;
    
    println!("\nGEV Parameters:");
    println!("  Location (μ): {:.3}", gev_params.location);
    println!("  Scale (σ): {:.3}", gev_params.scale);
    println!("  Shape (ξ): {:.3}", gev_params.shape);
    
    // Classify distribution type
    let dist_type = if gev_params.shape.abs() < 0.05 {
        "Gumbel (Type I)"
    } else if gev_params.shape > 0.0 {
        "Fréchet (Type II)"
    } else {
        "Weibull (Type III)"
    };
    println!("  Distribution type: {}", dist_type);
    
    // Calculate return values
    let return_periods = vec![10, 25, 50, 100, 200, 500, 1000];
    let return_values = gev_analysis.return_values(&return_periods, 1000);
    
    println!("\nReturn Values:");
    for rp in &return_periods {
        if let Some(ci) = return_values.get(rp) {
            println!("  {}-year: {:.2}m (90% CI: {:.2}-{:.2}m)",
                rp, ci.estimate, ci.lower, ci.upper);
        }
    }
    
    // Example 2: Peak-Over-Threshold Method
    println!("\n=== Peak-Over-Threshold Method ===");
    let mut pot_analysis = PotAnalysis::new(wave_data.clone());
    
    // Select threshold
    let threshold_selection = pot_analysis.select_threshold();
    println!("Recommended threshold: {:.2}m", threshold_selection.recommended);
    
    // Fit GPD
    let gpd_params = pot_analysis.fit_gpd()?;
    println!("\nGPD Parameters:");
    println!("  Threshold (u): {:.3}", gpd_params.threshold);
    println!("  Scale (σ): {:.3}", gpd_params.scale);
    println!("  Shape (ξ): {:.3}", gpd_params.shape);
    println!("  Rate (λ): {:.3} events/year", gpd_params.rate);
    
    // POT return values
    let pot_return_values = pot_analysis.return_values_pot(&return_periods);
    println!("\nPOT Return Values:");
    for rp in &return_periods {
        if let Some(&value) = pot_return_values.get(rp) {
            println!("  {}-year: {:.2}m", rp, value);
        }
    }
    
    // Compare methods
    println!("\n=== Method Comparison ===");
    println!("Return Period | Annual Maxima | POT");
    println!("--------------|---------------|-------");
    for rp in &[50, 100, 200] {
        let am_val = return_values.get(rp).map(|ci| ci.estimate).unwrap_or(0.0);
        let pot_val = pot_return_values.get(rp).unwrap_or(&0.0);
        println!("{:>13} | {:>13.2} | {:>6.2}", rp, am_val, pot_val);
    }
    
    // Risk assessment
    println!("\n=== Risk Assessment ===");
    let design_life = 50;  // years
    let design_levels = vec![
        ("Current 100-year design", return_values[&100].estimate),
        ("Climate-adjusted (+10%)", return_values[&100].estimate * 1.1),
        ("Conservative design", return_values[&200].estimate),
    ];
    
    for (name, level) in design_levels {
        let encounter_prob = calculate_encounter_probability(level, design_life, &gev_params);
        println!("{}: {:.2}m", name, level);
        println!("  Probability of exceedance in {} years: {:.1}%",
            design_life, encounter_prob * 100.0);
    }
    
    // Generate diagnostic plots
    eva_visualization::plot_threshold_diagnostics(
        &threshold_selection,
        "threshold_selection.png"
    )?;
    
    // Environmental considerations
    println!("\n=== Environmental Context ===");
    println!("- Climate change may invalidate stationarity assumption");
    println!("- Consider non-stationary EVA with time-varying parameters");
    println!("- Indigenous oral histories provide validation for extreme events");
    println!("- Compound events (waves + surge) require multivariate EVA");
    println!("- Precautionary principle suggests using upper confidence limits");
    
    Ok(())
}

/// Extract annual maxima from wave time series
fn extract_annual_maxima(data: &[(DateTime<Utc>, f64)]) -> Vec<f64> {
    use std::collections::HashMap;
    
    let mut annual_max: HashMap<i32, f64> = HashMap::new();
    
    for (timestamp, hs) in data {
        let year = timestamp.year();
        annual_max.entry(year)
            .and_modify(|max| *max = max.max(*hs))
            .or_insert(*hs);
    }
    
    let mut maxima: Vec<f64> = annual_max.values().cloned().collect();
    maxima.sort_by(|a, b| a.partial_cmp(b).unwrap());
    maxima
}

/// Calculate encounter probability
fn calculate_encounter_probability(
    level: f64,
    design_life: u32,
    params: &GevParams,
) -> f64 {
    let annual_prob = 1.0 - gev_cdf(level, params);
    1.0 - (1.0 - annual_prob).powi(design_life as i32)
}

/// Load wave data (placeholder)
fn load_wave_data() -> Result<Vec<(DateTime<Utc>, f64)>, Box<dyn Error>> {
    // In practice, load from file or database
    // For now, generate synthetic extremes
    use rand::{thread_rng, Rng};
    use rand_distr::{Distribution, LogNormal};
    
    let mut rng = thread_rng();
    let mut data = Vec::new();
    
    // Generate 20 years of 3-hourly data
    let start = Utc.ymd(2000, 1, 1).and_hms(0, 0, 0);
    for hours in (0..(20 * 365 * 24)).step_by(3) {
        let timestamp = start + chrono::Duration::hours(hours);
        
        // Log-normal distribution for wave heights
        let dist = LogNormal::new(0.5, 0.8).unwrap();
        let hs = dist.sample(&mut rng).min(15.0);  // Cap at 15m
        
        data.push((timestamp, hs));
    }
    
    Ok(data)
}
```

### Exercises

#### Exercise 1: Non-stationary EVA
Implement time-varying parameters for climate change:

```rust
pub struct NonStationaryGev {
    pub location_trend: f64,  // Linear trend in μ
    pub scale_trend: f64,     // Trend in σ
    pub base_year: i32,
}

impl NonStationaryGev {
    pub fn fit_with_covariates(
        &self,
        data: &[(i32, f64)],  // (year, annual_max)
    ) -> NonStationaryParams {
        // TODO: Implement MLE with time-varying parameters
        // μ(t) = μ₀ + μ₁ * (t - t₀)
        // σ(t) = σ₀ * exp(σ₁ * (t - t₀))
    }
}
```

#### Exercise 2: Multivariate Extremes
Analyze joint extremes of waves and water levels:

```rust
pub fn bivariate_pot_analysis(
    wave_surge_data: &[(f64, f64)],  // (Hs, surge)
) -> BivariateExtremes {
    // TODO: Implement bivariate threshold selection
    // Fit conditional extremes model
    // Calculate joint return periods
}
```

#### Exercise 3: Regional Frequency Analysis
Pool data from multiple sites for better estimates:

```rust
pub struct RegionalAnalysis {
    pub sites: Vec<SiteData>,
    pub homogeneity_test: HomogeneityResult,
}

impl RegionalAnalysis {
    pub fn index_flood_method(&self) -> RegionalGrowthCurve {
        // TODO: Implement L-moment based regional analysis
        // Test for regional homogeneity
        // Derive regional growth curve
    }
}
```

#### Exercise 4: Bayesian EVA
Implement Bayesian extreme value analysis:

```rust
pub fn bayesian_gev(
    data: &[f64],
    prior: GevPrior,
    n_mcmc: usize,
) -> PosteriorSamples {
    // TODO: Implement MCMC for GEV
    // Use informative priors from regional information
    // Generate posterior predictive distributions
}
```

### Questions for Reflection

1. **Stationarity**: Under what conditions might the stationarity assumption fail? How would you detect and account for trends in extreme wave heights?

2. **Threshold Selection**: What are the trade-offs in POT threshold selection? Too high loses data, too low violates GPD assumptions - how do you balance?

3. **Uncertainty Communication**: How would you communicate uncertainty in 100-year wave heights to stakeholders? Consider different risk tolerances.

4. **Compound Events**: Why might analyzing wave extremes in isolation underestimate coastal flood risk? How do waves, surge, and tide interact?

5. **Traditional Knowledge**: How can Indigenous oral histories of extreme events validate or challenge statistical extrapolations?

### Additional Resources

- Coles (2001) "An Introduction to Statistical Modeling of Extreme Values"
- IPCC Special Report on Extremes (SREX)
- R packages: extRemes, ismev, evd
- Python: pyextremes, climextRemes
- Coastal Engineering Manual - Chapter on Statistics

### Next Tutorial Preview

In Tutorial 4.3, we'll explore seasonal wave patterns, investigating how wave climates vary throughout the year and the implications for coastal processes, ecology, and human activities.

Remember: Extreme value analysis extrapolates beyond our observations. Always consider physical limits, communicate uncertainty, and design for resilience rather than just meeting standards.
