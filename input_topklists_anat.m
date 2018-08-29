%% input topklists
load medianvois % median values computed from AAL parcellation
medians_ctrl = medianvois(control_idx, :);
medians_sla = medianvois(sla_idx, :);
medians_park = medianvois(park_idx, :);
[Y,rank_sla] = sort(abs(medians_sla),2,'descend');
[Y,rank_park] = sort(abs(medians_park),2,'descend');
[Y,rank_ctrl] = sort(abs(medians_ctrl),2,'descend');
save md_anatvoi_abs_median_ranks rank_sla rank_park rank_ctrl

