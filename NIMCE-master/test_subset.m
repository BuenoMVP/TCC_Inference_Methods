% Test script with a subset of the data
% This will help verify the algorithm works before running on the full dataset

clear; clc;
fprintf('=== NIMCE Test with Data Subset ===\n');

% Load full dataset
data_path = '/home/marco/projects/TCC_Inference_Methods/Database/input data/net3_expression_data.tsv';
fprintf('Loading full dataset...\n');
full_data = dlmread(data_path, '\t');
fprintf('Full dataset: %d samples x %d genes\n', size(full_data,1), size(full_data,2));

% Create a subset for testing (first 200 samples, first 100 genes)
fprintf('Creating test subset: 200 samples x 100 genes\n');
test_data = full_data(1:200, 1:100);

% Parameters
timepoint = 20;  % 200 samples / 20 timepoints = 10 samples
window_num = 3;
time_lag = 2;
h = 0.45;
lambda2 = 0.04;

fprintf('Test parameters:\n');
fprintf('  - Timepoints: %d\n', timepoint);
fprintf('  - Window number: %d\n', window_num);
fprintf('  - Time lag: %d\n', time_lag);
fprintf('  - Kernel width: %.2f\n', h);
fprintf('  - Lambda: %.3f\n', lambda2);

% Run test
fprintf('\nRunning NIMCE on test subset...\n');
tStart = tic;

try
    [G] = NIMCE_optimized(test_data, timepoint, window_num, time_lag, h, lambda2);
    
    elapsed_time = toc(tStart);
    fprintf('\nTest completed successfully!\n');
    fprintf('Execution time: %.2f seconds (%.2f minutes)\n', elapsed_time, elapsed_time/60);
    fprintf('Result matrix: %d x %d\n', size(G,1), size(G,2));
    
    % Save test results
    test_output = '/home/marco/projects/TCC_Inference_Methods/NIMCE-master/NIMCE_test_subset.txt';
    dlmwrite(test_output, G, 'delimiter', '\t', 'precision', 100, 'newline', 'pc');
    fprintf('Test results saved to: %s\n', test_output);
    
catch ME
    fprintf('\nTest failed with error:\n');
    fprintf('Error: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
end

