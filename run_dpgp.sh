# conda deactivate
# conda activate py27

DP_GP_cluster.py \
-i "/nas/homes/benyang/Maresins/Tables/all_fc_zscore.tsv" \
--plot --plot_types pdf \
--output "/nas/homes/benyang/Maresins/DPGP_within_timepoint/all_within_timepoint"

# DP_GP_cluster.py \
# -i "/nas/homes/benyang/Maresins/Tables/sig_fc_zscore.tsv" \
# --plot --plot_types pdf \
# --output "/nas/homes/benyang/Maresins/DPGP_within_timepoint/sig_within_timepoint"