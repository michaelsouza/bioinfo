function x = matrix2array(x)
[nrow,ncol] = size(x);
if ncol > 1
    x = reshape(x, [nrow * 3,1]);
end
end