function x = array2matrix(x)
[nrow,ncol] = size(x);
if ncol == 1
    % convert array to matrix
    x = reshape(x, [nrow / 3, 3]);
end
end