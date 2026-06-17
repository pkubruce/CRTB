# CRTB

[中文](README.zh.md) | English

A Shiny-based interactive platform incorporating the extended China Tuberculosis Recurrence (CRTB) model.

This repository contains the code and data used for the extended China Tuberculosis Recurrence (CRTB) model analysis. The project is managed with the `renv` package (version 1.0.11). The complete R environment, including all required packages, is pre-packaged in the `renv/library/` folder — simply extract the archive and you are ready to run the code without any additional installation.

## Data

- All TB burden data were obtained from the World Health Organization (WHO) Global Tuberculosis Report 2025 and are stored in the `GTB2025/` folder.

## Code structure

The analysis pipeline is organized as follows:

- **`R/`** – Contains core utility functions and tools used throughout the analysis.
- **`01_fitted.R`** – Calibrates the model to WHO epidemiological targets.
- **`02_simulation.R`** – Runs baseline simulations for TB burden projections.
- **`03_Intervention.R`** – Simulates the four intervention scenarios (vaccination for high/low-risk recurrence groups, TPT for LTBI, and combined intervention).
- **`04_Manuscript.R`** – Produces the figures, tables, and summary statistics presented in the manuscript.

## Interactive Shiny application

The model functionality is also integrated into an interactive Shiny application. To run the app locally:

```r
shiny::runApp("app.R")
