import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

sample = ""
base_dir = ""
myRCTD_path = "cellbin_myRCTD.Rdata"
norm_weights_path = "optimized_norm_weights_mixed.rds"
print(f"Sample: {sample}")
print(f"Base directory: {base_dir}")

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()
numpy2ri.activate()

base = importr('base')
utils = importr('utils')
Matrix = importr('Matrix')

print("Reading norm_weights file...")
ro.r(f'norm_weights_raw <- readRDS("{norm_weights_path}")')
ro.r('norm_weights_mat <- as.matrix(norm_weights_raw)')

norm_weights_mat = ro.r('norm_weights_mat')
norm_weights_array = np.array(norm_weights_mat)
print(f"norm_weights array shape: {norm_weights_array.shape}")

print("Loading RCTD data...")
ro.r(f'load("{myRCTD_path}")')
ro.r('results <- myRCTD@results')
ro.r('spatialRNA <- myRCTD@spatialRNA')
ro.r('cell_type_names <- myRCTD@cell_type_info$info[[2]]')

coords_r = ro.r('spatialRNA@coords')
if hasattr(coords_r, 'shape'):
    coords = np.array(coords_r)
else:
    coords = np.array(coords_r)

cell_type_names_r = ro.r('cell_type_names')
cell_type_names = list(cell_type_names_r)
print(f"Original cell type names: {cell_type_names}")

print("Creating weights dataframe...")
n_columns = norm_weights_array.shape[1]
n_cell_types = len(cell_type_names)

if n_columns == n_cell_types + 1:
    print("Detected extra column, likely 'Mixed'. Adding to cell type names.")
    cell_type_names_with_mixed = cell_type_names + ["Mixed"]
    columns_to_use = cell_type_names_with_mixed
elif n_columns == n_cell_types:
    columns_to_use = cell_type_names
else:
    print(f"Warning: Column count mismatch! Using generic column names.")
    columns_to_use = [f"CellType_{i}" for i in range(n_columns)]

weights_df = pd.DataFrame(norm_weights_array, columns=columns_to_use)

if "Mixed" in weights_df.columns:
    weights_without_mixed = weights_df.drop("Mixed", axis=1)
    main_cell_types = weights_without_mixed.idxmax(axis=1)
    main_cell_type_probs = weights_without_mixed.max(axis=1)
    mixed_mask = weights_df["Mixed"] > 0.5
    if mixed_mask.any():
        print(f"Found {mixed_mask.sum()} spots with high Mixed probability")
        main_cell_types[mixed_mask] = "Mixed"
        main_cell_type_probs[mixed_mask] = weights_df.loc[mixed_mask, "Mixed"]
else:
    main_cell_types = weights_df.idxmax(axis=1)
    main_cell_type_probs = weights_df.max(axis=1)

print(f"Main cell types distribution:\n{main_cell_types.value_counts()}")

n_spots = norm_weights_array.shape[0]
if coords.shape[0] == 2 and coords.shape[1] == n_spots:
    coords = coords.T
elif coords.shape[0] == n_spots and coords.shape[1] == 2:
    pass
else:
    if coords.size == n_spots * 2:
        coords = coords.reshape(n_spots, 2)
    else:
        print(f"ERROR: Cannot reshape coords. coords.size={coords.size}, expected {n_spots*2}")

if coords.shape[0] == n_spots:
    spots_df = pd.DataFrame({
        'x': coords[:, 0],
        'y': coords[:, 1],
        'cell_type': main_cell_types.values,
        'cell_type_prob': main_cell_type_probs.values
    })
    print(f"Spots dataframe shape: {spots_df.shape}")
    print(f"\nCell type counts:\n{spots_df['cell_type'].value_counts()}")
else:
    min_rows = min(n_spots, coords.shape[0])
    coords = coords[:min_rows, :]
    main_cell_types = main_cell_types.iloc[:min_rows]
    main_cell_type_probs = main_cell_type_probs.iloc[:min_rows]
    spots_df = pd.DataFrame({
        'x': coords[:, 0],
        'y': coords[:, 1],
        'cell_type': main_cell_types.values,
        'cell_type_prob': main_cell_type_probs.values
    })

print("\n" + "="*60)
print("Using phenoptr to calculate distances to Epithelial cells...")
print("="*60)

print("Loading phenoptr package...")
try:
    phenoptr = importr('phenoptr')
    print("phenoptr package loaded successfully.")
except Exception as e:
    print(f"phenoptr not found: {e}")
    raise

csv_path = f"{base_dir}/spots_data_{sample}.csv"
spots_df.to_csv(csv_path, index=False)
print(f"Spots data saved to CSV: {csv_path}")

r_code = f'''
library(phenoptr)
library(dplyr)

spots_data <- read.csv("{csv_path}")

cat("Number of spots loaded:", nrow(spots_data), "\\n")
cat("Columns:", colnames(spots_data), "\\n")
cat("Unique cell types:", unique(spots_data$cell_type), "\\n")

if ("Mixed" %in% spots_data$cell_type) {{
  cat("Excluding Mixed cells from analysis...\\n")
  spots_data <- spots_data[spots_data$cell_type != "Mixed", ]
  cat("Number of spots after excluding Mixed:", nrow(spots_data), "\\n")
}}

if (!("Epithelial" %in% spots_data$cell_type)) {{
  cat("ERROR: No Epithelial cells found!\\n")
  quit()
}}

csd <- spots_data
colnames(csd)[colnames(csd)=="x"] <- "Cell X Position"
colnames(csd)[colnames(csd)=="y"] <- "Cell Y Position"
colnames(csd)[colnames(csd)=="cell_type"] <- "Phenotype"

csd$`Cell ID` <- seq_len(nrow(csd))
csd$`Field ID` <- "{sample}"

csd <- csd[, c(
  "Cell ID",
  "Field ID",
  "Cell X Position",
  "Cell Y Position",
  "Phenotype",
  "cell_type_prob"
)]

cat("\\nCell seg data prepared for phenoptr:\\n")
print(head(csd))

cat("\\nPhenotype counts in csd:\\n")
print(table(csd$Phenotype))

cat("\\nComputing nearest neighbor distances using phenoptr::find_nearest_distance...\\n")
distances_all <- find_nearest_distance(csd)
cat("\\n>>> Verifying phenoptr output...\\n")
cat(">>> Head of distances_all:\\n")
print(head(distances_all))
cat(">>> Column names of distances_all:\\n")
print(colnames(distances_all))
cat(">>> End of verification\\n\\n")
csd_with_dist <- bind_cols(csd, distances_all)

cat("\\nColumns after adding distances:\\n")
print(colnames(csd_with_dist))

epi_cols <- grep(
  "Epithelial",
  colnames(csd_with_dist),
  value = TRUE
)

if (length(epi_cols) == 0) {{
  stop("No Epithelial distance column found")
}}

cat("Using distance column:", epi_cols[1], "\\n")

csd_with_dist$distance_to_epi <- csd_with_dist[[epi_cols[1]]]

epi_cells <- csd_with_dist[csd_with_dist$Phenotype == "Epithelial", ]
other_cells <- csd_with_dist[csd_with_dist$Phenotype != "Epithelial", ]

cat("\\nNumber of Epithelial cells:", nrow(epi_cells), "\\n")
cat("Number of other cells:", nrow(other_cells), "\\n")

distance_summary <- other_cells %>%
  group_by(Phenotype) %>%
  summarise(
    mean_distance = mean(distance_to_epi, na.rm=TRUE),
    median_distance = median(distance_to_epi, na.rm=TRUE),
    sd_distance = sd(distance_to_epi, na.rm=TRUE),
    min_distance = min(distance_to_epi, na.rm=TRUE),
    max_distance = max(distance_to_epi, na.rm=TRUE),
    cell_count = n()
  ) %>%
  arrange(mean_distance)

cat("\\nDistance to nearest Epithelial cell by cell type:\\n")
print(distance_summary)

write.csv(distance_summary, file="{base_dir}/distance_to_epithelial_summary_{sample}.csv", row.names=FALSE)

other_cells_output <- other_cells[, c("Phenotype", "Cell X Position", "Cell Y Position",
                                      "cell_type_prob", "distance_to_epi")]
colnames(other_cells_output)[1] <- "cell_type"
colnames(other_cells_output)[2] <- "x"
colnames(other_cells_output)[3] <- "y"
write.csv(other_cells_output, file="{base_dir}/cell_distances_to_epithelial_{sample}.csv", row.names=FALSE)

list(
  distance_summary = distance_summary,
  cell_distances = other_cells_output,
  csd_with_dist = csd_with_dist
)
'''

try:
    print("Executing R code with phenoptr::find_nearest_distance...")
    result = ro.r(r_code)
    if result is not None and result is not ro.rinterface.NULL:
        distance_summary_r = result.rx2('distance_summary')
        distance_summary = pandas2ri.rpy2py(distance_summary_r)
        print(f"\nDistance summary DataFrame shape: {distance_summary.shape}")
        print("\nDistance to nearest Epithelial cell by cell type:")
        print(distance_summary.to_string())
        cell_distances_r = result.rx2('cell_distances')
        cell_distances = pandas2ri.rpy2py(cell_distances_r)
        print(f"\nDetailed distances DataFrame shape: {cell_distances.shape}")
        print(f"First few rows of cell distances:")
        print(cell_distances.head())
        print("\nVisualizing distance distributions...")
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        ax1 = axes[0, 0]
        cell_types_ordered = distance_summary['Phenotype']
        mean_distances = distance_summary['mean_distance']
        bars = ax1.bar(range(len(cell_types_ordered)), mean_distances)
        ax1.set_xlabel('Cell Type')
        ax1.set_ylabel('Mean Distance to Epi (pixels)')
        ax1.set_title(f'Mean Distance to Nearest Epithelial Cell - {sample}')
        ax1.set_xticks(range(len(cell_types_ordered)))
        ax1.set_xticklabels(cell_types_ordered, rotation=45, ha='right')
        for i, (cell_type, bar) in enumerate(zip(cell_types_ordered, bars)):
            if cell_type == 'Fib2':
                bar.set_color('red')
                ax1.text(i, bar.get_height() + 0.5, 'Fib2', ha='center', color='red', fontweight='bold')
        ax2 = axes[0, 1]
        ax2.hist(cell_distances['distance_to_epi'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        ax2.set_xlabel('Distance to Nearest Epi Cell (pixels)')
        ax2.set_ylabel('Frequency')
        ax2.set_title(f'Distribution of Distances to Epithelial Cells - {sample}')
        mean_all = cell_distances['distance_to_epi'].mean()
        ax2.axvline(mean_all, color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_all:.2f}')
        ax2.legend()
        ax3 = axes[1, 0]
        if 'Fib2' in cell_distances['cell_type'].unique():
            fib2_distances = cell_distances[cell_distances['cell_type'] == 'Fib2']['distance_to_epi']
            if len(fib2_distances) > 0:
                ax3.hist(fib2_distances, bins=30, alpha=0.7, color='lightcoral', edgecolor='black')
                ax3.set_xlabel('Distance to Nearest Epi Cell (pixels)')
                ax3.set_ylabel('Frequency')
                ax3.set_title(f'Fib2 Cells Distance to Epithelial - {sample}')
                fib2_mean = fib2_distances.mean()
                fib2_median = fib2_distances.median()
                ax3.axvline(fib2_mean, color='red', linestyle='--', linewidth=2,
                           label=f'Mean: {fib2_mean:.2f}')
                ax3.axvline(fib2_median, color='blue', linestyle=':', linewidth=2,
                           label=f'Median: {fib2_median:.2f}')
                ax3.legend()
            else:
                ax3.text(0.5, 0.5, 'No Fib2 cells found', ha='center', va='center', transform=ax3.transAxes)
                ax3.set_title('Fib2 Cells Distance to Epithelial (No Data)')
        else:
            ax3.text(0.5, 0.5, 'No Fib2 cells found', ha='center', va='center', transform=ax3.transAxes)
            ax3.set_title('Fib2 Cells Distance to Epithelial (No Data)')
        ax4 = axes[1, 1]
        scatter = ax4.scatter(cell_distances['distance_to_epi'],
                             cell_distances['cell_type_prob'],
                             alpha=0.5, s=10)
        ax4.set_xlabel('Distance to Nearest Epi Cell (pixels)')
        ax4.set_ylabel('Cell Type Probability')
        ax4.set_title(f'Distance vs Cell Type Probability - {sample}')
        if len(cell_distances) > 1:
            z = np.polyfit(cell_distances['distance_to_epi'],
                          cell_distances['cell_type_prob'], 1)
            p = np.poly1d(z)
            ax4.plot(cell_distances['distance_to_epi'],
                    p(cell_distances['distance_to_epi']),
                    color='red', linewidth=2,
                    label=f'Trend: y={z[0]:.4f}x+{z[1]:.4f}')
            ax4.legend()
        plt.tight_layout()
        plt.savefig(f"{base_dir}/distance_to_epithelial_analysis_{sample}.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{base_dir}/distance_to_epithelial_analysis_{sample}.pdf", bbox_inches='tight')
        plt.show()
        print("\nStatistical analysis of distances...")
        if 'Fib2' in cell_distances['cell_type'].unique():
            fib2_distances = cell_distances[cell_distances['cell_type'] == 'Fib2']['distance_to_epi']
            other_distances = cell_distances[cell_distances['cell_type'] != 'Fib2']['distance_to_epi']
            print(f"Fib2 cells: {len(fib2_distances)} cells, mean distance: {fib2_distances.mean():.2f}")
            print(f"Other cells: {len(other_distances)} cells, mean distance: {other_distances.mean():.2f}")
            from scipy.stats import mannwhitneyu
            u_stat, p_value = mannwhitneyu(
                fib2_distances,
                other_distances,
                alternative="two-sided"
            )
            print(f"\nMann-Whitney U test comparing Fib2 vs other cells:")
            print(f"  U-statistic: {u_stat:.4f}")
            print(f"  p-value: {p_value:.4f}")
            if p_value < 0.05:
                if fib2_distances.mean() < other_distances.mean():
                    print(f"  Conclusion: Fib2 cells are significantly closer to Epithelial cells than other cell types (p < 0.05)")
                else:
                    print(f"  Conclusion: Fib2 cells are significantly farther from Epithelial cells than other cell types (p < 0.05)")
            else:
                print(f"  Conclusion: No significant difference in distance to Epithelial cells between Fib2 and other cell types")
        print("\nSaving results...")
        distance_summary.to_csv(f"{base_dir}/distance_to_epithelial_summary_{sample}.csv", index=False)
        print(f"Distance summary saved to: {base_dir}/distance_to_epithelial_summary_{sample}.csv")
        cell_distances.to_csv(f"{base_dir}/cell_distances_to_epithelial_{sample}.csv", index=False)
        print(f"Detailed distance data saved to: {base_dir}/cell_distances_to_epithelial_{sample}.csv")
        with open(f"{base_dir}/distance_analysis_summary_{sample}.txt", "w") as f:
            f.write(f"Distance Analysis Summary - {sample}\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total spots analyzed: {len(spots_df)}\n")
            f.write(f"Epithelial cells: {(spots_df['cell_type'] == 'Epithelial').sum()}\n")
            if "Mixed" in spots_df['cell_type'].unique():
                f.write(f"Mixed cells: {(spots_df['cell_type'] == 'Mixed').sum()}\n")
            f.write(f"Non-epithelial cells (excluding Mixed): {(spots_df['cell_type'] != 'Epithelial').sum()}\n\n")
            f.write("Distance to nearest Epithelial cell by cell type:\n")
            f.write(distance_summary.to_string() + "\n\n")
            if 'Fib2' in cell_distances['cell_type'].unique():
                f.write("\nFib2-specific analysis:\n")
                f.write(f"  Number of Fib2 cells: {len(fib2_distances)}\n")
                f.write(f"  Mean distance to Epi: {fib2_distances.mean():.2f}\n")
                f.write(f"  Median distance to Epi: {fib2_distances.median():.2f}\n")
                f.write(f"  Min distance to Epi: {fib2_distances.min():.2f}\n")
                f.write(f"  Max distance to Epi: {fib2_distances.max():.2f}\n")
                f.write(f"  Mann-Whitney U p-value vs other cells: {p_value:.4f}\n")
                if p_value < 0.05:
                    if fib2_distances.mean() < other_distances.mean():
                        f.write("  Conclusion: Fib2 cells are significantly closer to Epithelial cells\n")
                    else:
                        f.write("  Conclusion: Fib2 cells are significantly farther from Epithelial cells\n")
                else:
                    f.write("  Conclusion: No significant difference in distance\n")
        print(f"Analysis summary saved to: {base_dir}/distance_analysis_summary_{sample}.txt")
    else:
        print("Distance calculation returned no results. Check if there are Epithelial cells.")
except Exception as e:
    print(f"Error in distance calculation: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*50)
print("DISTANCE ANALYSIS COMPLETED")
print("="*50)
print(f"Analysis saved in: {base_dir}")
print("="*50)