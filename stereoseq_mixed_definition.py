import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import stereo as st
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import cKDTree
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
import tempfile
import shutil

warnings.filterwarnings("ignore")
pandas2ri.activate()
numpy2ri.activate()

output_dir = ""
os.makedirs(output_dir, exist_ok=True)

samples = [
    {
        "name": "A06061C3_TN",
        "gef_path": "A06061C3.tissue.gef",
        "myRCTD_path": "A06061C3_TN_cellbin_myRCTD.Rdata",
        "norm_weights_path": "A06061C3_TN_cellbin.norm_weights.rds"
    },
    {
        "name": "A06061D1_TN",
        "gef_path": "A06061D1.tissue.gef",
        "myRCTD_path": "A06061D1_TN_cellbin_myRCTD.Rdata",
        "norm_weights_path": "A06061D1_TN_cellbin.norm_weights.rds"
    },
    {
        "name": "A06061G1_TN",
        "gef_path": "A06061G1.tissue.gef",
        "myRCTD_path": "A06061G1_TN_cellbin_myRCTD.Rdata",
        "norm_weights_path": "A06061G1_TN_cellbin.norm_weights.rds"
    },
    {
        "name": "A06265E5_TN",
        "gef_path": "A06265E5.tissue.gef",
        "myRCTD_path": "A06265E5_TN_cellbin_myRCTD.Rdata",
        "norm_weights_path": "A06265E5_TN_cellbin.norm_weights.rds"
    },
    {
        "name": "A06270A3_TN",
        "gef_path": "A06270A3.tissue.gef",
        "myRCTD_path": "A06270A3_TN_cellbin_myRCTD.Rdata",
        "norm_weights_path": "A06270A3_TN_cellbin.norm_weights.rds"
    },
    {
        "name": "A06270D1_TN",
        "gef_path": "A06270D1.tissue.gef",
        "myRCTD_path": "A06270D1_TN_cellbin_myRCTD.Rdata",
        "norm_weights_path": "A06270D1_TN_cellbin.norm_weights.rds"
    },
    {
        "name": "A06270E5_TN",
        "gef_path": "A06270E5.tissue.gef",
        "myRCTD_path": "A06270E5_TN_cellbin_myRCTD.Rdata",
        "norm_weights_path": "A06270E5_TN_cellbin.norm_weights.rds"
    },
    {
        "name": "Y01416J6_TN",
        "gef_path": "Y01416J6.tissue.gef",
        "myRCTD_path": "Y01416J6_TN_cellbin_myRCTD.Rdata",
        "norm_weights_path": "Y01416J6_TN_cellbin.norm_weights.rds"
    },
    {
        "name": "Y01416L6_TN",
        "gef_path": "Y01416L6.tissue.gef",
        "myRCTD_path": "Y01416L6_TN_cellbin_myRCTD.Rdata",
        "norm_weights_path": "Y01416L6_TN_cellbin.norm_weights.rds"
    },
    {
        "name": "Y01416M4_TN",
        "gef_path": "Y01416M4.tissue.gef",
        "myRCTD_path": "Y01416M4_TN_cellbin_myRCTD.Rdata",
        "norm_weights_path": "Y01416M4_TN_cellbin.norm_weights.rds"
    }
    ]

def run_stereopy_analysis(gef_path, bin_size=50):
    data = st.io.read_gef(gef_path, bin_size=bin_size, gene_name_index=True)
    data.tl.cal_qc()
    data.tl.filter_cells(min_counts=200, min_genes=200, pct_counts_mt=20)
    data.tl.normalize_total(target_sum=1e4)
    data.tl.log1p()
    data.tl.highly_variable_genes(n_top_genes=2000)
    data.tl.scale(zero_center=False)
    data.tl.pca(n_pcs=30)
    data.tl.neighbors(n_pcs=30)
    data.tl.spatial_neighbors()
    data.tl.leiden(neighbors_res_key="spatial_neighbors", res_key="spatial_leiden")
    return data, data.tl.result["spatial_leiden"]["group"].astype(str), data.position

def load_rctd_results(myRCTD_path, norm_weights_path):
    ro.r(f'load("{myRCTD_path}")')
    ro.r(f'norm_weights <- readRDS("{norm_weights_path}")')
    coords = np.array(ro.r("myRCTD@spatialRNA@coords"))
    weights = np.array(ro.r("as.matrix(norm_weights)"))
    cell_types = list(ro.r("myRCTD@cell_type_info$info[[2]]"))
    weights_df = pd.DataFrame(weights, columns=cell_types)
    main_type = weights_df.idxmax(axis=1)
    confidence = weights_df.max(axis=1).values
    return coords, weights_df, main_type, confidence

def align_coordinates(st_coords, r_coords, tol=5):
    tree = cKDTree(st_coords)
    d, idx = tree.query(r_coords, k=1)
    m = d <= tol
    return pd.DataFrame({"r": np.where(m)[0], "s": idx[m]})

def mark_mixed_by_structure(df):
    df = df.copy()
    df["final_type"] = df["cell_type"].copy()

    for cl in df["spatial_cluster"].unique():
        sub = df[df["spatial_cluster"] == cl]
        if len(sub) < 80:
            continue

        type_freq = sub["cell_type"].value_counts(normalize=True)
        dominant_type = type_freq.index[0]
        dominant_ratio = type_freq.iloc[0]

        if dominant_ratio > 0.85:
            continue

        coords = sub[["x", "y"]].values
        tree = cKDTree(coords)
        types = sub["cell_type"].values

        mixed_idx = []

        for i in range(len(sub)):
            _, nn = tree.query(coords[i:i+1], k=min(30, len(sub)))
            nn_types = types[nn[0]]
            nn_freq = pd.Series(nn_types).value_counts(normalize=True)

            local_dom = nn_freq.index[0]
            local_ratio = nn_freq.iloc[0]

            if (
                local_ratio < 0.65 and
                local_dom != dominant_type and
                sub.iloc[i]["confidence"] < 0.80
            ):
                mixed_idx.append(sub.index[i])

        df.loc[mixed_idx, "final_type"] = "Mixed"

    return df

def visualize(df, sample):
    color_mapping = {
        "Fib1": "#F5D2A8",
        "Fib2": "#FCED82", 
        "T_NK_cell": "#EE934E",
        "B_cell": "#3C77AF",
        "Epithelial": "#D1352B",
        "Plasma_cell": "#F5CFE4",
        "Macrophage": "#9B5B33",
        "cDCs": "#B383B9",
        "pDCs": "#8FA4AE",
        "EC": "#AECDE1",
        "LEC": "#D2EBC8",
        "Adipocyte": "#7F7F7F",
        "Mixed": "#CCCCCC"
    }
    
    plt.figure(figsize=(12, 10))
    
    all_types = df["final_type"].unique()
    
    ordered_types = []
    for cell_type in color_mapping.keys():
        if cell_type in all_types:
            ordered_types.append(cell_type)
    
    for cell_type in all_types:
        if cell_type not in color_mapping and cell_type not in ordered_types:
            ordered_types.append(cell_type)
    
    type_to_color = {}
    color_list = []
    type_list = []
    
    for cell_type in ordered_types:
        if cell_type in color_mapping:
            type_to_color[cell_type] = color_mapping[cell_type]
        else:
            type_to_color[cell_type] = "#999999"
        color_list.append(type_to_color[cell_type])
        type_list.append(cell_type)
    
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(color_list)
    
    type_to_idx = {cell_type: i for i, cell_type in enumerate(type_list)}
    color_indices = df["final_type"].map(type_to_idx).values
    
    scatter = plt.scatter(df["x"], df["y"], c=color_indices, cmap=cmap, s=6)
    
    plt.axis("equal")
    plt.title(f"{sample} - Mixed Filtered")
    
    handles = [
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=type_to_color[cell_type],
                   markersize=8)
        for cell_type in type_list
    ]
    
    plt.legend(handles, type_list,
               bbox_to_anchor=(1.02, 1),
               loc="upper left",
               frameon=False,
               title="Cell Types")
    
    plt.savefig(os.path.join(output_dir, f"{sample}_Mixed_filtered.pdf"),
                bbox_inches="tight")
    plt.close()

def create_rctd_probability_plots(sample_name, myRCTD_path, norm_weights_path):
    temp_dir = tempfile.mkdtemp()
    rctd_plots_dir = os.path.join(output_dir, f"{sample_name}_RCTD_plots")
    os.makedirs(rctd_plots_dir, exist_ok=True)
    
    try:
        r_code = f"""
        library(spacexr)
        library(Matrix)
        library(Seurat)
        library(ggplot2)
        library(RColorBrewer)

        setwd("{temp_dir}")
        load("{myRCTD_path}")
        norm_weights <- readRDS("{norm_weights_path}")
        cell_type_names <- myRCTD@cell_type_info$info[[2]]
        spatialRNA <- myRCTD@spatialRNA
        resultsdir <- "{temp_dir}"
        plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
        plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
        print("Generated files:")
        print(list.files(resultsdir, pattern="\\.pdf$"))
        """

        ro.r(r_code)

        for filename in os.listdir(temp_dir):
            if filename.endswith('.pdf'):
                source_path = os.path.join(temp_dir, filename)
                dest_path = os.path.join(rctd_plots_dir, filename)
                shutil.copy2(source_path, dest_path)
                print(f"Copied {filename} to {rctd_plots_dir}")

        self_contained_pdf_path = os.path.join(output_dir, f"{sample_name}_RCTD_probability_plots.pdf")
        self_contained_combined_pdf_path = os.path.join(output_dir, f"{sample_name}_RCTD_probability_plots_combined.pdf")

        self_contained_pdfs = []
        for filename in os.listdir(temp_dir):
            if filename.endswith('.pdf') and 'weights' in filename.lower():
                self_contained_pdfs.append(os.path.join(temp_dir, filename))

        if self_contained_pdfs:
            from PyPDF2 import PdfMerger
            merger = PdfMerger()
            
            for pdf_file in sorted(self_contained_pdfs):
                merger.append(pdf_file)
            
            merger.write(self_contained_combined_pdf_path)
            merger.close()
            
            print(f"Created combined PDF: {self_contained_combined_pdf_path}")

        prob_cmap = LinearSegmentedColormap.from_list(
            "probability_cmap",
            ["#0000FF", "#00BFFF", "#2E8B57", "#FFD700", "#FF0000"]
        )

        self_contained_all_types_pdf_path = os.path.join(output_dir, f"{sample_name}_cell_type_probabilities.pdf")

        norm_weights = ro.r(f'readRDS("{norm_weights_path}")')
        weights_matrix = np.array(ro.r("as.matrix(norm_weights)"))
        cell_types = list(ro.r("colnames(norm_weights)"))

        ro.r(f'load("{myRCTD_path}")')
        coords = np.array(ro.r("myRCTD@spatialRNA@coords"))

        with PdfPages(self_contained_all_types_pdf_path) as pdf:
            for i, cell_type in enumerate(cell_types):
                fig, ax = plt.subplots(figsize=(12, 10))

                weights = weights_matrix[:, i]

                scatter = ax.scatter(coords[:, 0], coords[:, 1], c=weights, 
                                     cmap=prob_cmap, s=6, vmin=0, vmax=1)
                
                ax.axis("equal")
                ax.set_title(f"{sample_name} - {cell_type} Probability Distribution")

                cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
                cbar.set_label('Prediction Probability', rotation=270, labelpad=15)

                mean_prob = np.mean(weights[weights > 0]) if np.any(weights > 0) else 0
                max_prob = np.max(weights) if len(weights) > 0 else 0
                cells_with_prob = np.sum(weights > 0)
                
                stats_text = f"Mean probability: {mean_prob:.3f}\nMax probability: {max_prob:.3f}\nCells with probability > 0: {cells_with_prob}"

                fig.text(0.02, 0.02, stats_text, transform=fig.transFigure,
                        verticalalignment='bottom', fontsize=10,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

                plt.subplots_adjust(bottom=0.15)

                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
        
        print(f"Created cell type probability plots for {sample_name}: {self_contained_all_types_pdf_path}")
        

        self_contained_threshold_pdf_path = os.path.join(output_dir, f"{sample_name}_cell_type_probabilities_threshold_0.25.pdf")

        with PdfPages(self_contained_threshold_pdf_path) as pdf:
            for i, cell_type in enumerate(cell_types):
                fig, ax = plt.subplots(figsize=(12, 10))

                weights = weights_matrix[:, i]

                threshold = 0.25
                threshold_mask = weights > threshold
                
                if np.any(threshold_mask):
                    scatter = ax.scatter(
                        coords[threshold_mask, 0], 
                        coords[threshold_mask, 1], 
                        c=weights[threshold_mask], 
                        cmap=prob_cmap, 
                        s=8,
                        vmin=0,
                        vmax=1
                    )
                    
                    ax.axis("equal")
                    ax.set_title(f"{sample_name} - {cell_type} Probability Distribution (Threshold > {threshold})")

                    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
                    cbar.set_label(f'Prediction Probability (>0.25)', rotation=270, labelpad=15)

                    threshold_weights = weights[threshold_mask]
                    mean_prob = np.mean(threshold_weights) if len(threshold_weights) > 0 else 0
                    max_prob = np.max(threshold_weights) if len(threshold_weights) > 0 else 0
                    cells_above_threshold = np.sum(threshold_mask)
                    
                    stats_text = (
                        f"Mean probability: {mean_prob:.3f}\n"
                        f"Max probability: {max_prob:.3f}\n"
                        f"Cells with probability > {threshold}: {cells_above_threshold} ({cells_above_threshold/len(weights)*100:.1f}%)\n"
                        f"Threshold applied: >{threshold}"
                    )
                else:
                    ax.axis("equal")
                    ax.set_title(f"{sample_name} - {cell_type} Probability Distribution (Threshold > {threshold})")
                    ax.text(0.5, 0.5, f"No cells with probability > {threshold}", 
                            horizontalalignment='center', verticalalignment='center',
                            transform=ax.transAxes, fontsize=12)
                    
                    stats_text = f"No cells with probability > {threshold}"

                fig.text(0.02, 0.02, stats_text, transform=fig.transFigure,
                        verticalalignment='bottom', fontsize=10,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

                plt.subplots_adjust(bottom=0.15)

                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
        
        print(f"Created thresholded cell type probability plots for {sample_name}: {self_contained_threshold_pdf_path}")
        
    except Exception as e:
        print(f"Error creating RCTD probability plots for {sample_name}: {str(e)}")
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def process_sample(s):
    st_data, spatial_cluster, st_coords = run_stereopy_analysis(s["gef_path"])
    r_coords, weights_df, main_type, confidence = load_rctd_results(
        s["myRCTD_path"], s["norm_weights_path"])

    match = align_coordinates(st_coords, r_coords)
    df = pd.DataFrame({
        "x": st_coords[match.s, 0],
        "y": st_coords[match.s, 1],
        "spatial_cluster": spatial_cluster.iloc[match.s].values,
        "cell_type": main_type.iloc[match.r].values,
        "confidence": confidence[match.r]
    })

    df2 = mark_mixed_by_structure(df)
    visualize(df2, s["name"])

    opt = weights_df.copy()
    mixed_r = match.loc[df2["final_type"] == "Mixed", "r"]

    opt.iloc[mixed_r] = 0
    if "Mixed" not in opt.columns:
        opt["Mixed"] = 0
    opt.loc[mixed_r, "Mixed"] = 1

    ro.globalenv["opt"] = pandas2ri.py2rpy(opt)
    ro.r(f'saveRDS(opt, "{output_dir}/{s["name"]}_optimized_norm_weights_mixed.rds")')
    df2.to_csv(f"{output_dir}/{s['name']}_final_annotation.csv", index=False)

    n_total = len(df2)
    n_mixed = (df2["final_type"] == "Mixed").sum()
    print(f"{s['name']}: Total cells = {n_total}, Mixed cells = {n_mixed} ({n_mixed/n_total*100:.1f}%)")

    print(f"Creating RCTD probability plots for {s['name']}...")
    create_rctd_probability_plots(s["name"], s["myRCTD_path"], s["norm_weights_path"])

def main():
    for s in samples:
        print(f"Processing {s['name']}...")
        process_sample(s)
        print(f"Completed {s['name']}\n")

if __name__ == "__main__":
    main()