function [G] = NIMCE_optimized(data, timepoint, window_num, time_lag, h, lambda2)
% NIMCE Optimized with progress tracking
% Handles large datasets more efficiently

fprintf('Step 1/4: Preparing data...\n');
[row, ~] = size(data);
nsample = row / timepoint;

fprintf('  - Total samples: %d\n', nsample);
fprintf('  - Timepoints per sample: %d\n', timepoint);
fprintf('  - Creating time-lagged windows...\n');

window_data = createTime_lag(data, nsample, window_num-1);
fprintf('  - Time-lagged windows created successfully\n');

fprintf('\nStep 2/4: Computing zero-order interactions...\n');
fprintf('  - This is the most computationally intensive step\n');
fprintf('  - Progress will be shown for each window combination\n');

G_zero_order = [];
for ws = 1:window_num
    for we = 1:window_num
        if (we - ws) == time_lag
            fprintf('  - Processing window combination: ws=%d, we=%d\n', ws, we);
            
            % Show progress for this window combination
            G_zero_order = creatG_optimized(window_data{ws}, window_data{ws+1}, ...
                                          window_data{we}, window_data{ws+1}, ...
                                          window_data{we}, h);
            fprintf('  - Zero-order matrix computed: %dx%d\n', size(G_zero_order,1), size(G_zero_order,2));
        end
    end
end

fprintf('\nStep 3/4: Edge reduction...\n');
fprintf('  - Applying threshold filtering...\n');

G = edgereduce_optimized(data, G_zero_order, timepoint, window_num, time_lag, h, lambda2);

fprintf('\nStep 4/4: Final processing...\n');
fprintf('  - Computing final network...\n');

% Final network computation
nsample = row / timepoint;
window_data = createTime_lag(data, nsample, window_num-1);

for ws = 1:window_num
    for we = 1:window_num
        if (we - ws) == time_lag
            fprintf('  - Final network computation: ws=%d, we=%d\n', ws, we);
            G = creatGcmi(G, window_data{ws}, window_data{ws+1}, window_data{ws+2}, h, lambda2);
        end
    end
end

fprintf('  - Network inference completed!\n');
end

