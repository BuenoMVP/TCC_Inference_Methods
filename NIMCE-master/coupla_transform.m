function [ x_index_i ] = coupla_transform(x)
% Copula-transform array - rank and scale to [0, 1]
    [~, x_index] = sort(x, 2);
    [~, x_index_i] = sort(x_index, 2);
    x_index_i = (x_index_i - 1) / (size(x_index_i, 2) - 1);
end

