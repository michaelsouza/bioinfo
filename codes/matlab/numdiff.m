function g = numdiff(f, x)
fx = f(x);
nfun = length(fx);
d = zeros(size(x));
h = 0.0001;
g = zeros(nfun, length(x));
for k = 1:length(x)
    d(k) = h * max(abs(x(k)), 1.0);
    fx = f(x + d) - f(x - d);
    g(:,k) = fx / (2 * d(k));
    d(k) = 0.0;
end
if nfun == 1
    % gradient
    g = g';
end
end

function diff()
end

function grad()
end

function jacobian()
end