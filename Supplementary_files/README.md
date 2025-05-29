# Metabolite Comparison using Venn Diagrams and Bar Plots

This complementary seccion contains an R script for visualizing shared and unique precursor ion features across experimental groups from metabolomics data. The analysis was performed for both positive and negative ionization in both adquisitions modes of ions (5-15) using Venn diagrams and bar plots.

## Description

The code accomplishes the following:

1. **Reads Excel data** with multiple sheets (positive and negative ionization modes).
2. **Filters data by experimental group** based on ionization energy and method.
3. **Generates three pairwise Venn diagrams** for each ionization mode:
   - 5_m vs 5_sm
   - 15_m vs 15_sm
   - 5_m vs 15_m
4. **Calculates shared features (intersection counts)** for each comparison.
5. **Creates a bar plot** showing the total and shared precursor counts.
6. **Combines all figures** (Venn diagrams and bar plots) into a single composite figure using `gridExtra` and `ggpubr`.
7. **Applies color schemes**:
   - Shades of red for positive mode
   - Shades of blue for negative mode
   - Matched ionization subgroups share hues with differing intensity
   - Cross-comparison (5_m vs 15_m) is highlighted with a neutral color

## Final Figure

The final figure integrates all visualizations into a single multi-panel layout. Each subplot is labeled (A–H) and grouped by ionization mode.

(Result/venn_results/Figuras/PNG/figure_union.png))

## Dependencies
- R (>= 4.0.0)
- Packages:
  - `readxl`
  - `tidyverse`
  - `VennDiagram`
  - `ggplot2`
  - `gridExtra`
  - `ggpubr`
  
## Suggested Usage
This script can be adapted for:
- Metabolomics data exploration
- Quality control in feature detection workflows
- Pre-processing diagnostics for multigroup comparisons
