function problem = create_mdgp_instance_rand(n, dij_max, dij_eps, omega_vals, omega_prob)

if nargin < 2 || isempty(dij_max)
    dij_max = 5;
end

if nargin < 3 || isempty(dij_eps)
    dij_eps = 0;
end

if nargin < 4 || isempty(omega_vals)
    omega_vals = [60, 180, 300];
    omega_prob = [1/3,1/3,1/3];
end

% convert to radians
omega_vals = pi * omega_vals / 180;
omega_eps = 5 * pi / 180; % 5 grades


% bond lengths
d = 1.526 * ones(n,1);

% set planar angles
theta = 1.91 * ones(n,1);

% set dihedral angles
omega_prob = cumsum(omega_prob);
omega = zeros(n,1);
for k = 1:n
    omega(k) = omega_vals(sum(rand > omega_prob) + 1);
end
omega = omega + (2 * rand(n,1) - 1) * omega_eps;

% calculates euclidean coordinates
x = calc_x_from_w(d,theta,omega);

% create edges
i = zeros(n^2,1);
j = zeros(n^2,1);
dij_eps = dij_eps * ones(n^2,1);
nedges = 0;
for xi = 1:n
    for xj = (xi+1):n
        if xj - 1 <= 2
            nedges = nedges + 1;
            dij_eps(nedges) = 0; 
            i(nedges) = xi;
            j(nedges) = xj;
        elseif norm(x(xi,:) - x(xj,:)) < dij_max
            nedges = nedges + 1;
            i(nedges) = xi;
            j(nedges) = xj;
        end
    end
end
i = i(1:nedges);
j = j(1:nedges);

problem = create_mdgp_instance(x,i,j,dij_eps);
end