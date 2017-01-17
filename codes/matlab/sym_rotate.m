clear z y p theta x

%% create symbolic variables
z = sym('z', [3, 1], 'real');
y = sym('y', [3, 1], 'real');
p = sym('p', [3, 1], 'real');
theta = sym('theta');

%% unpacking
d = (y - z);     % direction
d = d / norm(d); % normalizing
y1=y(1); y2=y(2); y3=y(3); % base point
d1=d(1); d2=d(2); d3=d(3); % line direction
p1=p(1); p2=p(2); p3=p(3); % point to be rotated

%% evaluate simple expressions
sint = sin(theta);
cost = cos(theta);
omcost = 1 - cost;
sum_pd = sum(p.*d);

%% create auxiliary vectors
C = [(y1*(d2^2+d3^2) - d1*(y2*d2+y3*d3-sum_pd)); ...\
     (y2*(d1^2+d3^2) - d2*(y1*d1+y3*d3-sum_pd)); ...
     (y3*(d1^2+d2^2) - d3*(y1*d1+y2*d2-sum_pd))];
S = [(-y3*d2 + y2*d3 - d3*p2 + d2*p3);...
    ( y3*d1 - y1*d3 + d3*p1 - d1*p3);...
    (-y2*d1 + y1*d2 - d2*p1 + d1*p2)];

%% rotation
x  = C * omcost + p .* cost + S * sint;

%% diffs
g.theta = diff(x,theta);
g.p = jacobian(x,p);
g.y = jacobian(x,y);
g.z = jacobian(x,z);

%% create functions
fprintf('Creating rotate.m file\n')
y = [x,g.theta,g.p,g.y,g.z];
rotate = matlabFunction(y,'File','rotate.m');