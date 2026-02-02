# ExtremeTracker.jl

[![Julia](https://img.shields.io/badge/Julia-1.6+-9558B2?style=flat&logo=julia&logoColor=white)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive Julia package for tracking and analyzing marine/atmospheric heat wave events in spatiotemporal data. This package implements three state-of-the-art algorithms for heat wave detection, tracking, and composite analysis.

## Overview

**HeatWaveTracker.jl** provides three complementary methods for heat wave analysis:

1. **`hwtrack_nouniform`** - Spatially Coherent Tracking (Sun et al., 2023)
2. **`Tracker`** (Ocetrac) - Spatiotemporally Coherent Tracking (Scannell et al., 2023)
3. **`SpatialTemporalNormalization`** - Spatial-temporal normalization for composite analysis (Zhao et al., in review)

These methods can be used independently or combined to provide comprehensive heat wave event characterization from detection through composite analysis.

## Features

- 🔍 **Multiple Tracking Algorithms**: Choose from two state-of-the-art tracking methods
- 🌊 **Handles Complex Events**: Track splitting, merging, and evolution of heat wave events
- 📊 **Composite Analysis**: Normalize events to standard circular grids for inter-comparison
- 🚀 **High Performance**: Written in Julia for computational efficiency
- 🔧 **Flexible**: Works with various gridded temperature datasets (SST, 2m temperature, etc.)
- 📈 **Rich Metrics**: Extract intensity, duration, area, and trajectory information

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/ZijieZhaoMMHW/HeatWaveTracker.jl")
```

Or in the Julia REPL:

```julia
] add https://github.com/ZijieZhaoMMHW/HeatWaveTracker.jl
```

### Dependencies

```julia
using Images
using ImageMorphology
using ImageSegmentation
using Statistics
using Interpolations
using LinearAlgebra
```

## Quick Start

### Method 1: Event Tracking with Splitting/Merging (Sun et al., 2023)

This method tracks heat waves while explicitly handling event splitting and merging.

```julia
using HeatWaveTracker

# Load your temperature anomaly data (3D: lon, lat, time)
# hw: binary mask of heat wave conditions (0/1)
# lon_grid: longitude coordinates
# lat_grid: latitude coordinates

# Track heat waves
tracks = hwtrack_nouniform(hw, lon_grid, lat_grid, nums=10)

# Each track contains:
# - day: time indices
# - xloc: x-coordinates (cell arrays for each time step)
# - yloc: y-coordinates (cell arrays for each time step)
# - ori_day: origin day
# - ori_order: origin order
# - split_num: number of splits
# - split_day: days when splits occurred
```

### Method 2: Object-Based Tracking (Scannell et al., 2023)

This method uses morphological operations and 3D connectivity for robust tracking.

```julia
using HeatWaveTracker

# Create tracker object
tracker = Tracker(
    data,                    # 3D temperature array (time, y, x)
    mask,                    # Binary mask (1=valid ocean, 0=land)
    radius=8,                # Morphological radius
    min_size_quartile=0.75,  # Size threshold (75th percentile)
    positive=true            # Track positive anomalies
)

# Run tracking
labels, metadata = track(tracker)

# Results:
# - labels: 3D array with unique IDs for each heat wave
# - metadata: dictionary with tracking statistics
println("Tracked $(metadata["final_objects_tracked"]) heat wave events")
```

### Method 3: Spatial-Temporal Normalization (Zhao et al., in review)

Normalize tracked events to a standard circular grid for composite analysis.

```julia
using HeatWaveTracker

# Step 1: Calculate maximum radius for each track
radius_max = calculate_max_radius(tracks, lon_grid, lat_grid)

# Step 2: Normalize to circular grid
normalized_data = spatial_temporal_normalization(
    tracks,           # Tracked events from Method 1 or 2
    data_anomaly,     # 4D anomaly data (lon, lat, time, variables)
    lon_grid,         # Longitude coordinates
    lat_grid,         # Latitude coordinates
    radius_max,       # Maximum radii
    res=50,           # Spatial resolution (50x50 grid)
    n_phases=5        # Temporal phases (0, 0.25, 0.5, 0.75, 1.0)
)

# Step 3: Compute composite mean
composite = compute_composite(normalized_data, method="mean")

# Output dimensions: (res, res, n_phases, n_variables)
```


## Method Comparison

| Feature | Sun et al. (2023) | Scannell et al. (2023) |
|---------|-------------------|------------------------|
| **Splitting/Merging** | Explicit tracking | Implicit via connectivity |
| **Morphological Ops** | No | Yes (closing/opening) |
| **Size Filtering** | Fixed threshold | Percentile-based |
| **Best For** | Complex evolution | Clean boundaries |
| **Speed** | Slow | Fast |
| **Memory** | Moderate | Small |

**Recommendation**: 
- Use **Sun et al.** for tracking events with complex splitting/merging behavior
- Use **Scannell et al.** for cleaner event boundaries and detailed morphological filtering
- Both methods produce compatible output for normalization

## Key Parameters

### hwtrack_nouniform
- `nums`: Minimum number of pixels for valid heat wave (default: 10)
- `para_alpha`: Overlap threshold for tracking (default: 0.5)

### Tracker (Ocetrac)
- `radius`: Morphological structuring element radius (recommended: 4-10)
- `min_size_quartile`: Size threshold as quantile (0-1, recommended: 0.75)
- `positive`: Track positive (true) or negative (false) anomalies

### SpatialTemporalNormalization
- `res`: Spatial resolution of normalized grid (default: 50)
- `n_phases`: Number of temporal phases (default: 5)

## Output Data Structures

### Track Structure
```julia
struct Track
    day::Vector{Int}                    # Time indices
    xloc::Vector{Vector{Float32}}       # X-coordinates per timestep
    yloc::Vector{Vector{Float32}}       # Y-coordinates per timestep
    ori_day::Int                        # Origin day
    ori_order::Int                      # Origin order
    split_num::Vector{Int}              # Split counts
    split_day::Vector{Int}              # Split days
end
```

### Normalized Data
5D Array: `(spatial_x, spatial_y, temporal_phase, track_id, variable)`
- Dimensions 1-2: Normalized circular grid (polar coordinates)
- Dimension 3: Lifecycle phase (0=genesis, 1=termination)
- Dimension 4: Individual track/event ID
- Dimension 5: Variable (e.g., SST, SLP, wind)

## Applications

This package has been used for:
- 🌊 Marine heat wave tracking and characterization
- 🌡️ Atmospheric heat wave analysis  
- 🌀 Composite analysis of extreme events
- 📊 Statistical analysis of heat wave properties
- 🔬 Climate model validation

## Performance Tips

1. **Use appropriate data resolution**: 0.25° to 1° works well for global analysis
2. **Adjust parameters based on scale**: Larger radius for larger events
3. **Pre-filter small events**: Use `nums` or `min_size_quartile` to reduce memory
4. **Parallel processing**: Julia's built-in parallelism can speed up tracking
5. **Memory management**: Process time chunks for very long datasets

## Citation

If you use this package, please cite the relevant methods:

**For hwtrack_nouniform:**
```bibtex
@article{sun2023,
  title={Characterizing global marine heatwaves under a spatio-temporal framework},
  author={Sun, D. and Jing, Z. and Li, F. and Wu, L.},
  journal={Progress in Oceanography},
  volume={211},
  pages={102947},
  year={2023}
}
```

**For Tracker (Ocetrac):**
```bibtex
@article{scannell2023,
  title={Spatiotemporal Evolution of Marine Heatwaves Globally},
  author={Scannell, H. A. and Cai, C. and Thompson, L. and Whitt, D. B. and Gagne, D. J. and Abernathey, R.},
  journal={Journal of Geophysical Research: Oceans},
  year={2023}
}
```




## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes:

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Development

### Running Tests
```julia
using Pkg
Pkg.test("HeatWaveTracker")
```

### Building Documentation
```julia
using Pkg
Pkg.activate("docs")
include("docs/make.jl")
```

## Troubleshooting

**Common Issues:**

1. **Memory errors with large datasets**: Process data in chunks
2. **No events detected**: Check threshold and `min_size_quartile` parameters
3. **Tracking stops prematurely**: Verify data continuity and mask consistency
4. **Interpolation errors in normalization**: Ensure sufficient valid data points

See [Issues](https://github.com/yourusername/HeatWaveTracker.jl/issues) for more help.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.



## Related Projects

- [marineHeatWaves (Python)](https://github.com/ecjoliver/marineHeatWaves)
- [heatwaveR (R)](https://github.com/robwschlegel/heatwaveR)
- [Ocetrac (Python)](https://github.com/ocetrac/ocetrac)

---

**Keywords**: heat waves, marine heat waves, extreme events, climate, tracking, composite analysis, spatiotemporal analysis, Julia

