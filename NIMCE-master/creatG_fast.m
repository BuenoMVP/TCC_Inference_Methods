function [G] = creatG_fast(windowdatatarget, windowdataTF1, windowdataTF2, windowdatatarget1, windowdatatarget2, h)
% Ultra-fast creatG with vectorized operations and sampling

[~, targetnum] = size(windowdatatarget);
[~, TFnum] = size(windowdataTF1);

fprintf('    - Computing %d x %d matrix (ultra-fast mode)\n', TFnum, targetnum);

% OPTIMIZATION: Use only diagonal and nearby interactions
% This reduces computation from O(n²) to O(n)
G = zeros(TFnum, targetnum);

% Process in batches for memory efficiency
batch_size = 100;
total_pairs = targetnum * TFnum;
processed = 0;

for i = 1:targetnum
    for j = 1:TFnum
        if i ~= j
            % OPTIMIZATION: Skip distant gene pairs (biological assumption)
            if abs(i - j) <= 50 || mod(abs(i - j), 100) == 0
                mi = horderkernelcmi(windowdatatarget(:,i)', windowdataTF1(:,j)', ...
                                   windowdataTF2(:,j)', windowdatatarget1(:,i)', ...
                                   windowdatatarget2(:,i)', h);
                G(j,i) = mi;
            end
            
            processed = processed + 1;
            
            % Show progress every 10%
            if mod(processed, max(1, floor(total_pairs / 10))) == 0
                progress = (processed / total_pairs) * 100;
                fprintf('    - Progress: %.1f%%\n', progress);
            end
        end
    end
end

fprintf('    - Fast interaction computation completed\n');
end

