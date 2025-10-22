function [G] = edgereduce_fast(data, G_zero_order, timepoint, window_num, time_lag, h, lambda2)
% Ultra-fast edgereduce with optimized operations

[row, col] = size(data);

fprintf('    - Fast threshold computation...\n');

% OPTIMIZATION: Use quantile-based threshold instead of k-means
G_1d = sort(G_zero_order(:));
threshold_percentile = 0.95;  % Keep top 5% of connections
threshold_idx = round(length(G_1d) * threshold_percentile);
lamda1 = G_1d(threshold_idx);

fprintf('    - Threshold (95th percentile): %.6f\n', lamda1);

% Fast vectorized thresholding
G_zero_order(G_zero_order < lamda1) = 0;

% Count edges
non_zero_edges = sum(G_zero_order(:) > 0);
total_edges = numel(G_zero_order);
fprintf('    - Edges retained: %d/%d (%.2f%%)\n', non_zero_edges, total_edges, (non_zero_edges/total_edges)*100);

% Fast final processing
nsample = row / timepoint;
window_data = createTime_lag(data, nsample, window_num-1);

fprintf('    - Fast final network computation...\n');
for ws = 1:window_num
    for we = 1:window_num
        if (we - ws) == time_lag
            G = creatGcmi(G_zero_order, window_data{ws}, window_data{ws+1}, window_data{ws+2}, h, lambda2);
        end
    end
end

fprintf('    - Fast edge reduction completed\n');
end

