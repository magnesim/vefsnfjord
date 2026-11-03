# CAMABUR - Coastal Marine Radionuclide Bioaccumulation Model

Analysis and modeling tools for coastal marine transport and biological uptake of radionuclides (specifically $^{137}\text{Cs}$). Developed and validated in the Vefsnfjorden case study (Nordland, Norway). Full documentation is available here (link will come). 

---

## Overview

The repository couples physical ocean transport simulations with a radioecological bioaccumulation model. It simulates the transfer of radionuclides from seawater to marine biota at designated monitoring stations, accounting for hydrodynamics, dynamic kinetics, and salinity variations.

### Key Modeling Features

1. **Biological Uptake Schemes**:
   - **Dynamic Kinetic Model**: Solves first-order transfer kinetics:
     $$\frac{dC_b}{dt} = k_{\text{up}} C_w(t) - k_{\text{el}} C_b(t)$$
     where $k_{\text{el}} = \frac{\ln(2)}{t_{1/2}}$, and uptake rate $k_{\text{up}} = C_F \cdot k_{\text{el}}$.
   - **Instantaneous Model**: Equilibrium approximation $C_b(t) = C_F \cdot C_w(t)$.
   - **Salinity-Dependent Concentration Factors ($C_F(S)$)**: Stepwise empirical $C_F$ and $k_{\text{up}}$ scaled by ambient seawater salinity (e.g., from ROMS ocean simulations).

2. **Seawater Concentration Sources ($C_w$)**:
   - **OpenDrift**: Extracts concentration time series from OpenDrift dispersion netCDF outputs (`concentration_smooth`), supporting dissolved and total fractions.
   - **Trajan**: Reads trajectory/gridded concentration files via xarray.
   - **Synthetic/Analytical**: Generates baseline or idealized time series (`constant`, `decay`, `sinus`, `randomspikes`, `setzero`) for sensitivity tests.

3. **Spatial Mapping & Statistics**:
   - Matches observation coordinates to nearest ocean/transport model grid cells within a user-defined impact radius.
   - Correlates seawater salinity and radionuclide concentrations across stations.

---

## Repository Structure

```text
vefsnfjord/
├── bio_model/
│   ├── bio_uptake02.py        # Main execution driver
│   ├── bio_uptake_tools.py    # Uptake solvers, data extraction, and plotting
│   └── bio_input.yaml         # Configuration file
├── LICENSE
└── README.md
```

---

## Installation

We recommend using [Mamba](https://mamba.readthedocs.io/) or [Conda](https://docs.conda.io/) to handle C-library geospatial dependencies (`cartopy`, `netCDF4`).

```bash
# Clone the repository
git clone https://github.com/magnesim/vefsnfjord.git
cd vefsnfjord

# Create environment from environment.yaml
conda env create -f environment.yaml

# Activate the environment
conda activate camabur
```

---

## Configuration (`bio_input.yaml`)

Edit `bio_model/bio_input.yaml` to set your run parameters:

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `tstart`, `tend` | `str` | Simulation start and end timestamps (`YYYY-MM-DD HH`) |
| `model_folder`, `model_fn` | `str` | Path and filename pattern for transport NetCDF (uses `EXPNM` template) |
| `exp_nm` | `str` | Experiment / run label (e.g., `run2022-001`) |
| `water_conc` | `str` | Source: `model`, `constant`, `decay`, `sinus`, `randomspikes`, `setzero` |
| `bio_uptake_model` | `str` | Bio model formulation: `dynamic` or `instant` |
| `concentration_factor` | `float` | Equilibrium concentration factor $C_F$ ($\text{Bq/kg fw}$ per $\text{Bq/L}$) |
| `bio_half_life` | `float` | Biological half-life $t_{1/2}$ in days |
| `initial_bio_conc` | `float` | Initial biota activity concentration ($\text{Bq/kg fw}$) |
| `background_value` | `float` | Seawater background offset ($\text{Bq/L}$) |
| `salt_CF` | `bool` | Enable salinity-dependent $C_F$ lookups |
| `salt_intervals_CF` | `list` | Intervals `[min_salt, max_salt, CF]` |
| `salt_datadir`, `salt_filenames` | `str` | Path/glob to ROMS ocean history netCDF files |
| `stations_to_eval` | `list` | Station identifiers to evaluate |
| `stations_coords` | `dict` | Coordinates `[latitude, longitude, sample_time]` for stations |
| `output_folder` | `str` | Directory to save generated figures and CSV data |
| `saveplots` | `bool` | Save output plots to disk (`True`) or display interactively (`False`) |

---

## Usage

Navigate to the `bio_model` directory and run `bio_uptake02.py`:

```bash
cd bio_model

# Run with default config (bio_input.yaml)
python bio_uptake02.py

# Run with a custom config file and verbose logging
python bio_uptake02.py -f /path/to/custom_config.yaml --verbose
```

### CLI Arguments

- `-f`, `--configfile`: Path to the YAML configuration file (default: `bio_input.yaml`).
- `--verbose`: Enable detailed processing and extraction logs.

---

## Outputs

When `saveplots: True` is configured:

1. **Tabular Results (`.csv`)**:
   - `bio_model_results_{exp_nm}_00.csv`: Matched observations containing station name, coordinates, sample timestamp, simulated seawater $^{137}\text{Cs}$ concentration ($\text{Bq/L}$), and biota $^{137}\text{Cs}$ activity concentration ($\text{Bq/kg fresh weight}$).
2. **Visualizations (`.png`)**:
   - `bio_uptake_map_{exp_nm}.png`: Station locations and bathymetry/topography contours.
   - `bio_model_{exp_nm}06_dissolved.png`: Time-series comparison of seawater concentration, biota concentration, apparent $C_F$, and salinity.
   - `bio_model_{exp_nm}_stats.png`: Histograms and scatter correlation between salinity and radionuclide concentrations.

---

## License

This project is licensed under the terms of the license file included in the repository (`LICENSE`).
