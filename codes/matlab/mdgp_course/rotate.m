function y = rotate(theta,a,b,x)
% Rotates x around the axis ab with angle theta
% theta: rotation angle
% a : axis start
% b : axis end
% x : point to be rotated

% converting to columns
x = x(:);
a = a(:);
b = b(:);

% get direction
d = (b - a);     % direction
d = d / norm(d); % normalizing

% using the notation of Murray's website
A=b(1); B=b(2); C=b(3); % base point
U=d(1); V=d(2); W=d(3); % line direction
X=x(1); Y=x(2); Z=x(3); % point to be rotated

% evaluate simple expressions
dot_xd = sum(x.*d);

% create auxiliary vectors
u = [(A*(V^2+W^2) - U*(B*V+C*W-dot_xd)); ...\
     (B*(U^2+W^2) - V*(A*U+C*W-dot_xd)); ...
     (C*(U^2+V^2) - W*(A*U+B*V-dot_xd))];
 
v = [(-C*V + B*W - W*Y + V*Z);...
    ( C*U - A*W + W*X - U*Z);...
    (-B*U + A*V - V*X + U*Y)];

% rotation
y  = u * (1 - cos(theta)) + x * cos(theta) + v * sin(theta);