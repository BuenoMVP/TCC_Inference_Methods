function [G] = edgereduce_optimized(data, G_zero_order, timepoint, window_num, time_lag, h, lambda2)
% Optimized version of edgereduce with progress tracking

[row, col] = size(data);

fprintf('    - Applying k-means clustering for threshold determination...\n');
G_1d = sort(G_zero_order(:));

% Use a subset for k-means if the matrix is too large
if length(G_1d) > 10000
    fprintf('    - Large matrix detected, using subset for k-means...\n');
    subset_size = min(10000, length(G_1d));
    subset_indices = round(linspace(1, length(G_1d), subset_size));
    G_subset = G_1d(subset_indices);
    [~, Ctrs] = kmeans(G_subset, 2, 'Start', 'uniform', 'Replicates', 3);
else
    [~, Ctrs] = kmeans(G_1d, 2, 'Start', 'uniform', 'Replicates', 3);
end

lamda1 = (Ctrs(1) + Ctrs(2)) / 2;
fprintf('    - Threshold lambda1: %.6f\n', lamda1);

fprintf('    - Filtering edges below threshold...\n');
% Vectorized thresholding for better performance
G_zero_order(G_zero_order < lamda1) = 0;

% Count non-zero edges
non_zero_edges = sum(G_zero_order(:) > 0);
total_edges = numel(G_zero_order);
fprintf('    - Edges retained: %d/%d (%.2f%%)\n', non_zero_edges, total_edges, (non_zero_edges/total_edges)*100);

nsample = row / timepoint;
window_data = createTime_lag(data, nsample, window_num-1);

fprintf('    - Computing final network with CMI...\n');
for ws = 1:window_num
    for we = 1:window_num
        if (we - ws) == time_lag
            fprintf('    - Final CMI computation: ws=%d, we=%d\n', ws, we);
            G = creatGcmi(G_zero_order, window_data{ws}, window_data{ws+1}, window_data{ws+2}, h, lambda2);
        end
    end
end

fprintf('    - Edge reduction completed\n');
end

