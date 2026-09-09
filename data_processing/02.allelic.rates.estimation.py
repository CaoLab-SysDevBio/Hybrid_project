# ============================================================
# we used three replicates at LatTII stage as an example.
# We choose scDALI "Joint" model to identify genes with allelic imbalance among different cell states. 
# ============================================================

### first step: identify genes with allelic imbalance ###
import os
os.chdir("/allelic_expression_matrices/")
print(os.getcwd()) #view current path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scdali import run_tests, run_interpolation
from scdali.utils import simulate
from scdali.utils.stats import apply_fdr_bh, compute_quantile_diff

replicates = ["LTB2_Cs1Ci1", "LTB2_Cs1Ci2", "LTB2_Cs2Ci1"] 

for rep in replicates:
    print(f"\nRunning replicate: {rep}")
    total_file = f"/{rep}_cutoff30.counts.f.total.txt"
    maternal_file = f"/{rep}_cutoff30.counts.f.maternal.txt"
    pca_file = f"exp_seurat_obj/{rep}_cutoff30_PCA_dim50_coordinates.csv"
    # input count of D and A - cell x gene
    counts_OG_all = pd.read_csv(total_file, sep="\t", index_col=0)
    print(type(counts_OG_all))
    counts_OG_all.shape
    counts_OG_maternal=pd.read_csv(maternal_file, sep="\t", index_col=0)
    print(type(counts_OG_maternal))
    counts_OG_maternal.shape
    # input cell state
    pc_cor = pd.read_csv(pca_file, sep=",", index_col=0)
    cell_state_pc = pc_cor.iloc[:, :20]
    # for CsCi RNA: genexcell
    counts_OG_maternal_A = counts_OG_maternal.T
    counts_OG_all_D = counts_OG_all.T
    # delete gene: A>D
    bad_OG_mask = (counts_OG_maternal_A > counts_OG_all_D).any(axis=0)
    bad_OGs = counts_OG_maternal_A.columns[bad_OG_mask]
    print(f"Number of bad OGs: {len(bad_OGs)}")
    # keep only valid OGs
    counts_OG_maternal_A = counts_OG_maternal_A.loc[:, ~bad_OG_mask]
    counts_OG_all_D = counts_OG_all_D.loc[:, ~bad_OG_mask]
    print(counts_OG_maternal_A.shape, counts_OG_all_D.shape, cell_state_pc.shape)
    ## joint - for each gene
    success_joint_pc = {}
    error_records_joint_pc = []
    for col in counts_OG_maternal_A.columns:
        try:
            res = run_tests(A=counts_OG_maternal_A[col], D=counts_OG_all_D[col], cell_state=cell_state_pc, model='scDALI-Joint',base_rate=0.5, return_rho=True, n_cores=10)
            success_joint_pc[col] = {"pv_joint": res['pvalues'][0], "rho": res['rhos'][0]}
        except Exception as e:
            error_records_joint_pc.append({
                "OG": col,
                "error_type": type(e).__name__,
                "error_message": str(e)
            })
            print(f"Joint failed: {col}")
            continue
    pv_rho_joint_df_pc = pd.DataFrame.from_dict(success_joint_pc, orient='index').reset_index()
    pv_rho_joint_df_pc.columns = ["OG", "pv_joint", "rho"]
    error_joint_pc_df = pd.DataFrame(error_records_joint_pc)
    pv_rho_joint_df_pc.to_csv(f"/{rep}_cutoff30_joint_pc20_pv_rho_all.csv", index=False)
    error_joint_pc_df.to_csv(f"/{rep}_cutoff30_joint_pc20_error_info.csv", index=False)
    print(f"Finished: {rep}")


### second step: calculate allelic-bias rates ###
import os
import sys
os.chdir("/allelic_expression_matrices/")
print(os.getcwd())

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scdali import run_tests, run_interpolation
from scdali.utils import simulate
from scdali.utils.stats import apply_fdr_bh, compute_quantile_diff

# if len(sys.argv) != 2:
#     raise ValueError("Usage: python run_interpolation.py <replicate_name>")
# rep = sys.argv[1]
# print("=" * 60)
# print(f"Running replicate: {rep}")
# print("=" * 60)

total_file = f"/{rep}_cutoff30.counts.f.total.txt"
maternal_file = f"/{rep}_cutoff30.counts.f.maternal.txt"
pca_file = f"exp_seurat_obj/{rep}_cutoff30_PCA_dim50_coordinates.csv"
# input count of D and A - cell x gene
counts_OG_all = pd.read_csv(total_file, sep="\t", index_col=0)
counts_OG_maternal=pd.read_csv(maternal_file, sep="\t", index_col=0)
# input cell state
pc_cor = pd.read_csv(pca_file, sep=",", index_col=0)
cell_state_pc = pc_cor.iloc[:, :20]
# for CsCi RNA: genexcell
counts_OG_maternal_A = counts_OG_maternal.T
counts_OG_all_D = counts_OG_all.T
# pc-only sig p of three models
pv_pc = pd.read_csv(f"allelic_rates/{rep}_cutoff30_joint_sig_padj_0.01_OGlist.csv")
significant_OG_pc = pv_pc["OG"]
results_pc = run_interpolation(
    A = counts_OG_maternal_A.loc[:, significant_OG_pc],
    D = counts_OG_all_D.loc[:, significant_OG_pc],
    cell_state = cell_state_pc,
    maxiter = 3000,
    n_cores = 20
)
# calculate posterior_mean
posterior_mean_pc = results_pc['posterior_mean']
posterior_var_pc = results_pc["posterior_var"]
rate_means_pc = pd.DataFrame(posterior_mean_pc, columns = significant_OG_pc, index = cell_state_pc.index)
rate_vars_pc = pd.DataFrame(posterior_var_pc, columns = significant_OG_pc, index = cell_state_pc.index)
rate_means_pc.to_csv(f"/{rep}_cutoff30_joint_posterior_mean_pv0.01.csv")
rate_vars_pc.to_csv(f"/{rep}_cutoff30_joint_posterior_vars_pv0.01.csv")
# print("=" * 60)
# print(f"Finished successfully: {rep}")
# print("=" * 60)
