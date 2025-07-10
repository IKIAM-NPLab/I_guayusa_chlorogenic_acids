# Metabolite Comparison using Bar Plots

This complementary seccion contains an R script for visualizing shared and unique precursor ion features across experimental groups from metabolomics data. The analysis was performed for both positive and negative ionization in both adquisitions modes of ions (5-15) using bar plots.

## Description

The code accomplishes the following:

1. **Reads Excel data** with multiple sheets (positive and negative ionization modes).
2. **Filters data by experimental group** based on ionization energy and method.
3. **Generates groups** for each ionization mode:
   - 5_m (FastDDA 5 ions with merge) vs 5_sm (FastDDA 5 ions without merge)
   - 15_m (FastDDA 15 ions with merge) vs 15_sm (FastDDA 15 ions without merge)
4. **Calculates shared features (intersection counts)** for each comparison.
5. **Creates a bar plot** showing the total and shared precursor counts.
6. **Applies color schemes**:
   - Shades of red for positive mode
   - Shades of blue for negative mode
   - Matched ionization subgroups share hues with differing intensity

![Final Composite Figure](Result/venn_results/Figuras/PNG/figure_union_2.png)

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
