% NIMCE Optimized Script with Progress Indicators
% Modified to handle large datasets with progress tracking

clear; clc;
fprintf('=== NIMCE Algorithm - Optimized Version ===\n');
fprintf('Loading dataset...\n');

% Timer start
tStart = tic;

% Dataset configuration
comparedata = {'net3_expression_data'};
data_path = '/home/marco/projects/TCC_Inference_Methods/Database/input data/net3_expression_data.tsv';

% Load TSV data with progress indicator
fprintf('Reading data from: %s\n', data_path);
data = dlmread(data_path, '\t');
fprintf('Data loaded: %d samples x %d genes\n', size(data,1), size(data,2));

% Parameters
window_num = 3;
time_lag = 2;
h = 0.45;
timepoint = 13;  % 806 samples / 13 timepoints = 62 samples
lambda2 = 0.04;

% Check if data dimensions are compatible
if mod(size(data,1), timepoint) ~= 0
    error('Number of samples must be divisible by timepoint');
end

fprintf('Parameters:\n');
fprintf('  - Timepoints: %d\n', timepoint);
fprintf('  - Window number: %d\n', window_num);
fprintf('  - Time lag: %d\n', time_lag);
fprintf('  - Kernel width: %.2f\n', h);
fprintf('  - Lambda: %.3f\n', lambda2);

% Run NIMCE with progress tracking
fprintf('\nStarting NIMCE algorithm...\n');
fprintf('This may take several hours for large datasets.\n');
fprintf('Progress will be shown for each major step.\n\n');

try
    [G] = NIMCE_optimized(data, timepoint, window_num, time_lag, h, lambda2);
    
    % Save results
    output_file = ['/home/marco/projects/TCC_Inference_Methods/NIMCE-master/database', 'NIMCE_', comparedata{1}, '.txt'];
    fprintf('Saving results to: %s\n', output_file);
    dlmwrite(output_file, G, 'delimiter', '\t', 'precision', 100, 'newline', 'pc');
    
    fprintf('\n=== Algorithm completed successfully! ===\n');
    fprintf('Results saved to: %s\n', output_file);
    fprintf('Matrix dimensions: %d x %d\n', size(G,1), size(G,2));
    
catch ME
    fprintf('\n=== ERROR OCCURRED ===\n');
    fprintf('Error message: %s\n', ME.message);
    fprintf('Error location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    rethrow(ME);
end

% Show execution time
elapsed_time = toc(tStart);
fprintf('\nTotal execution time: %.2f seconds (%.2f minutes)\n', elapsed_time, elapsed_time/60);
