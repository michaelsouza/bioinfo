function [value, Q, x, y] = rmsd(x, y, idx, stitle)
% idx: only x(idx,:) and y(idx,:) coordinates are used in the computation
% x  : real matrix with three columns
% y  : real matrix with three columns

% check input parameters
[~,ncols] = size(x);
if(ncols < 3)
    error('RMSD::The input x must be a matrix with 3 columns\n');
end
[~,ncols] = size(y);
if(ncols < 3)
    error('RMSD::The input y must be a matrix with 3 columns\n');
end

% extracting the coords specified by idx
if(nargin > 2 && ~isempty(idx))
    x = x(idx,:);
    y = y(idx,:);
end

% get centers
n  = length(ix);
xc = sum(x) / length(x); % center of x
yc = sum(y) / length(y); % center of y

% translation
for i = 1:n
    x(i,:) = x(i,:) - xc;
    y(i,:) = y(i,:) - yc;
end

% svd decomposition
A = x' * y;
[U, ~, V] = svd(A);

% rotation
Q = U * V';
x = x * Q;

% set rmsd value
value = norm(x - y, 'fro') / sqrt(n);

% plot superimposed structures
if(nargin > 3 && ~isempty(stitle))
    stitle = sprintf('%s, value: %e, npoints: %d', stitle, value, n);
    view_coords(x,y,stitle);
end
end

function view_coords(x,y,stitle)
n = length(x);

figure;
hold on;
for i =	1:n
    plot3( x(i,1), x(i,2), x(i,3), '*b');
    plot3( y(i,1), y(i,2), y(i,3), 'or');
end
hold off;

if(nargin == 3)
    title(stitle);
else
    title('view coords');
end
end