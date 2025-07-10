# Coastal Dynamics Tutorial Series: Rust Implementation
## Focus on Vancouver and Saint-Lawrence Coastal Systems

### Course Overview
This tutorial series implements key concepts from Bosboom & Stive's "Coastal Dynamics" using Rust, with applications relevant to MarineLabs' coastal intelligence platform. We'll build tools for wave analysis, coastal erosion modelling, and port operations optimization.

### Learning Objectives
By completing this series, you will:
1. **Understand** fundamental coastal processes including waves, currents, and sediment transport
2. **Implement** numerical models for coastal engineering applications in Rust
3. **Analyze** real coastal data from Canadian waters (Vancouver and Saint-Lawrence)
4. **Visualize** coastal phenomena using Rust plotting libraries
5. **Apply** knowledge to real-world problems in marine safety and coastal resilience
6. **Connect** engineering solutions to environmental and social impacts

### Tutorial Structure

**Important: All tutorials are designed to be independent and self-contained.** While they are numbered to suggest a logical learning progression, you can complete them in any order based on your interests and needs. Each tutorial includes all necessary code and data structures, so no previous tutorials are required.

#### Chapter 3: Ocean Waves (4 tutorials)
- **3.1** Wave Statistics and Rust Fundamentals (Week 1)
- **3.2** Spectral Analysis Implementation (Week 2)
- **3.3** Wave Generation and Dispersion (Week 3)
- **3.4** Tide Prediction Engine (Week 4)

#### Chapter 4: Global Wave Environments (3 tutorials)
- **4.1** Wave Climate Analysis for BC Coast (Week 5)
- **4.2** Extreme Value Statistics (Week 6)
- **4.3** Seasonal Wave Patterns (Week 7)

#### Chapter 5: Coastal Hydrodynamics (5 tutorials)
- **5.1** Wave Transformation Modelling (Week 8)
- **5.2** Wave-Induced Currents (Week 9)
- **5.3** Tidal Propagation in Channels (Week 10)
- **5.4** Combined Wave-Current Interactions (Week 11)
- **5.5** Real-time Hydrodynamic Monitoring (Week 12)

#### Chapter 7: Cross-shore Transport (3 tutorials)
- **7.1** Beach Profile Evolution (Week 13)
- **7.2** Storm Erosion Modelling (Week 14)
- **7.3** Dune Impact Assessment (Week 15)

#### Chapter 8: Longshore Transport (3 tutorials)
- **8.1** Longshore Transport Calculations (Week 16)
- **8.2** Coastline Evolution Modelling (Week 17)
- **8.3** Port Sedimentation Analysis (Week 18)

#### Chapter 9: Coastal Inlets and Tidal Basins (4 tutorials)
- **9.1** Inlet Stability Analysis (Week 19)
- **9.2** Tidal Prism Calculations (Week 20)
- **9.3** Ebb-Delta Dynamics (Week 21)
- **9.4** Basin Morphology Evolution (Week 22)

#### Chapter 10: Coastal Protection (4 tutorials)
- **10.1** Coastal Vulnerability Assessment (Week 23)
- **10.2** Nature-Based Solutions Design (Week 24)
- **10.3** Climate Adaptation Strategies (Week 25)
- **10.4** Integrated Coastal Management System (Week 26)

### Suggested Reading Schedule

#### Phase 1: Intensive Learning (July-August, 8 weeks)
**Target: 5-10 hours/week**

| Week | Chapter Reading | Tutorial | Hours |
|------|----------------|----------|--------|
| 1 | Ch. 3.1-3.4 | Tutorial 3.1 | 8 |
| 2 | Ch. 3.5-3.6 | Tutorial 3.2 | 8 |
| 3 | Ch. 3.7-3.8 | Tutorial 3.3 | 8 |
| 4 | Ch. 3.9 | Tutorial 3.4 | 8 |
| 5 | Ch. 4.1-4.3 | Tutorial 4.1 | 8 |
| 6 | Ch. 4.4 | Tutorial 4.2 | 8 |
| 7 | Review Ch. 3-4 | Tutorial 4.3 | 8 |
| 8 | Ch. 5.1-5.2 | Tutorial 5.1 | 8 |

#### Phase 2: Steady Progress (September-December, 18 weeks)
**Target: 5 hours/week**

| Week | Chapter Reading | Tutorial | Hours |
|------|----------------|----------|--------|
| 9-12 | Ch. 5.3-5.8 | Tutorials 5.2-5.5 | 20 |
| 13-15 | Ch. 7 | Tutorials 7.1-7.3 | 15 |
| 16-18 | Ch. 8 | Tutorials 8.1-8.3 | 15 |
| 19-22 | Ch. 9 | Tutorials 9.1-9.4 | 20 |
| 23-26 | Ch. 10 | Tutorials 10.1-10.4 | 20 |

### Development Environment Setup

#### Required Tools
```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Verify installation
rustc --version
cargo --version

# Create project structure
cargo new coastal_dynamics --bin
cd coastal_dynamics
```

#### Key Dependencies
Add to `Cargo.toml`:
```toml
[dependencies]
# Numerical computing
ndarray = "0.15"
nalgebra = "0.32"
num-complex = "0.4"

# Data processing
polars = { version = "0.35", features = ["lazy", "temporal", "strings"] }
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"

# Visualization
plotters = "0.3"
plotly = "0.18"

# Scientific computing
rustfft = "6.1"
statrs = "0.16"

# HTTP requests for data
reqwest = { version = "0.11", features = ["json"] }
tokio = { version = "1", features = ["full"] }

# Error handling
anyhow = "1.0"
thiserror = "1.0"
```

### Data Sources

#### Vancouver Coast
- **Environment Canada Marine Weather**: https://weather.gc.ca/marine/
- **DFO Water Levels**: https://www.tides.gc.ca/
- **Ocean Networks Canada**: https://www.oceannetworks.ca/

#### Saint-Lawrence
- **SLGO (St. Lawrence Global Observatory)**: https://ogsl.ca/
- **Canadian Hydrographic Service**: https://www.charts.gc.ca/

### Connection to MarineLabs Applications

Each tutorial will relate to MarineLabs' core services:
1. **Real-time monitoring**: Building streaming data processors
2. **Wave analysis**: Implementing spectral analysis and statistics
3. **Vessel wake detection**: Modelling wave propagation
4. **Coastal resilience**: Erosion and flood risk assessment
5. **Port operations**: Tidal predictions and current analysis

### Environmental and Social Context

Throughout the tutorials, we'll consider:
- **Climate change impacts** on coastal communities
- **Indigenous coastal management** practices
- **Ecosystem services** provided by coastal environments
- **Sustainable development** strategies
- **Environmental justice** in coastal planning

### Mathematical Notation Standards

All tutorials consistently use LaTeX notation for mathematical symbols and equations. This ensures:
- **Clarity**: Mathematical concepts are unambiguous
- **Consistency**: Same notation throughout all tutorials
- **Professional presentation**: Aligns with academic standards
- **Readability**: Complex equations are properly formatted

Examples:
- Wave dispersion relation: $\omega^2 = gk \tanh(kh)$
- Significant wave height: $H_s = 4\sqrt{m_0}$
- JONSWAP spectrum: $S(f) = \alpha g^2 (2\pi)^{-4} f^{-5} \exp\left(-\frac{5}{4}\left(\frac{f_p}{f}\right)^4\right) \cdot \gamma^r$

This notation standard makes the mathematical content accessible and professional.

### Assessment Strategy

Each tutorial includes:
1. **Conceptual Questions**: Test understanding of coastal processes
2. **Coding Exercises**: Progressive Rust implementations
3. **Data Analysis**: Work with real coastal datasets
4. **Visualization Tasks**: Create meaningful plots
5. **Integration Projects**: Optional projects that combine concepts across chapters (for those who complete multiple tutorials)

### Getting Started

You can begin with any tutorial that interests you, as each is self-contained. However, if you're new to both Rust and coastal dynamics, Tutorial 3.1 provides a good introduction to both. Before starting any tutorial, ensure you have:
1. Completed the development environment setup
2. The relevant textbook chapter available for reference
3. Allocated appropriate time (see schedule above)

Each tutorial will provide or generate its own datasets, so you don't need to worry about dependencies between tutorials.

Let's build tools that contribute to safer, more resilient coastlines!
