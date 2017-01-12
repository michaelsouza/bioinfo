function g = numdiff(f, x)
fx = f(x);
[nrows,ncols] = size(fx);
g = zeros(nrows,ncols,length(x));
d = zeros(size(x));
h = 0.001;
for k = 1:length(x)
    d(k) = h * max(abs(x(k)), 1.0);
    fx = f(x + d) - f(x - d);
    g(:,:,k) = fx / (2 * d(k));
    d(k) = 0.0;
end