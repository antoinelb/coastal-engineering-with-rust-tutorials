# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Architecture

This is a Rust-based coastal dynamics simulation toolkit implementing wave statistics and visualization for educational purposes. The project structure follows a modular design:

### Core Modules
- **`src/main.rs`**: Entry point that runs wave analysis scenarios (calm, moderate, storm conditions)
- **`src/wave_stats.rs`**: Core wave analysis functionality including:
  - `WaveTimeSeries`: Time series data structure with wave detection algorithms
  - `WaveStatistics`: Statistical calculations (Hs, Hmax, Tmean)
  - Zero-crossing wave detection and port safety assessment
- **`src/data_generator.rs`**: Synthetic wave data generation using superposition of sinusoids
- **`src/visualization.rs`**: Plotting capabilities using the `plotters` crate for time series and histograms

### Key Dependencies
- `plotters`: For generating wave plots and histograms
- `rand` + `rand_distr`: For stochastic wave generation with normal distribution noise
- Edition 2024 with clippy linting configured

### Data Flow
1. Generate synthetic wave data for different sea conditions
2. Detect individual waves using zero-crossing analysis
3. Calculate wave statistics (significant height, maximum height, mean period)
4. Assess port safety based on wave height thresholds
5. Generate visualizations (time series plots and height distributions)

The application outputs PNG plots for each scenario and displays statistical summaries including environmental context related to climate change impacts on coastal systems.

## Claude Code Guidelines

- Never generate code unless specifically asked for
- Act as a helper only, providing guidance and support
