function [G] = creatG_optimized(windowdatatarget, windowdataTF1, windowdataTF2, windowdatatarget1, windowdatatarget2, h)
% Optimized version of creatG with progress tracking

[~, targetnum] = size(windowdatatarget);
[~, TFnum] = size(windowdataTF1);

fprintf('    - Computing %d x %d interaction matrix\n', TFnum, targetnum);
fprintf('    - Total interactions to compute: %d\n', TFnum * targetnum);

G = zeros(TFnum, targetnum);
total_interactions = TFnum * targetnum;
processed = 0;

% Process in batches to show progress
batch_size = max(1, floor(total_interactions / 20)); % Show progress every 5%

for i = 1:targetnum
    for j = 1:TFnum
        if i ~= j
            % Compute mutual information
            mi = horderkernelcmi(windowdatatarget(:,i)', windowdataTF1(:,j)', ...
                               windowdataTF2(:,j)', windowdatatarget1(:,i)', ...
                               windowdatatarget2(:,i)', h);
            G(j,i) = mi;
        end
        
        processed = processed + 1;
        
        % Show progress every batch_size interactions
        if mod(processed, batch_size) == 0 || processed == total_interactions
            progress = (processed / total_interactions) * 100;
            fprintf('    - Progress: %.1f%% (%d/%d interactions)\n', progress, processed, total_interactions);
        end
    end
end

fprintf('    - Interaction matrix computation completed\n');
end

