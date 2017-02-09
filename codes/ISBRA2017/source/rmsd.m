function [value, Q, x, y] = rmsd(x, y, idx, stitle)
% idx: only x(idx,:) and y(idx,:) coordinates are used in the computation
% x  : real matrix with three columns
% y  : real matrix with three columns

% check input parameters
x = x';
y = y';
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
n  = length(x);
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
    stitle = sprintf('%s: RMSD=%3.2E, #Points: %d', stitle, value, n);
    view_coords(x,y,stitle,leg);
end
end

function view_coords(x,y,stitle)
fig = figure;
n = length(x);
subplot(1,2,1);
view([32 16]);
hold on;
box  on;
grid on;
% Create axes
d = zeros(n,1);
for i =	1:n
    plot3(x(i,1), x(i,2), x(i,3), '*b');
    plot3(y(i,1), y(i,2), y(i,3), 'or');
    d(i) = norm(x(i,:) - y(i,:));
    if i == 1
        legend('x','y');
    end
end
hold off;

if(nargin > 2)
    title(stitle);
else
    title('view coords');
end
subplot(1,2,2)
plot(sort(d),'-x');
xlabel('Pair ID')
ylabel('d(x_i,y_i)')
title('RMSD: Distance Between Correspondent Points');
fig.Position = [95,380,1420,420];
end