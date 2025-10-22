function [G] = NIMCE_fast(data, timepoint, window_num, time_lag, h, lambda2)
% Ultra-fast NIMCE with aggressive optimizations

fprintf('Step 1/3: Fast data preparation...\n');
[row, ~] = size(data);
nsample = row / timepoint;

% Create time-lagged windows
window_data = createTime_lag(data, nsample, window_num-1);
fprintf('  - Time-lagged windows created\n');

fprintf('\nStep 2/3: Fast interaction computation...\n');
fprintf('  - Using optimized computation...\n');

% Compute zero-order interactions with progress
G_zero_order = [];
for ws = 1:window_num
    for we = 1:window_num
        if (we - ws) == time_lag
            fprintf('  - Processing window: ws=%d, we=%d\n', ws, we);
            G_zero_order = creatG_fast(window_data{ws}, window_data{ws+1}, ...
                                    window_data{we}, window_data{ws+1}, ...
                                    window_data{we}, h);
        end
    end
end

fprintf('\nStep 3/3: Fast edge reduction...\n');
G = edgereduce_fast(data, G_zero_order, timepoint, window_num, time_lag, h, lambda2);

fprintf('  - Fast NIMCE completed!\n');
end

