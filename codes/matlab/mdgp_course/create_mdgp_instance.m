function problem = create_mdgp_instance(x, i, j, dij_eps)
if nargin < 2 || isempty(i) || isempty(j)
    % creates a complete instance with all edges
    natoms = size(x, 1); % number of atoms
    nedges = natoms*(natoms-1)/2;
    i = zeros(nedges,1);
    j = zeros(nedges,1);
    nedges = 0; % number of edges
    for xi = 1:natoms
        for xj = (xi+1):natoms
            nedges = nedges + 1;
            i(nedges) = xi;
            j(nedges) = xj;
        end
    end
end
nedges = length(i);

if nargin < 4 || isempty(dij_eps)
    % set distance precision
    dij_eps = zeros(nedges, 1);
end

% calculates distances and constraints
d = zeros(nedges, 1);
l = zeros(nedges, 1);
u = zeros(nedges, 1);
for k = 1:nedges
    xi = i(k);
    xj = j(k);
    d(k) = norm(x(xi,:) - x(xj,:));
    l(k) = d(k) * (1 - dij_eps(k));
    u(k) = d(k) * (1 + dij_eps(k));
end
problem.nedges = nedges;
problem.natoms = size(x,1);
problem.x = x; % euclidean coordinates
problem.i = i;
problem.j = j;
problem.d = d; % exact distance
problem.l = l; % lower bound
problem.u = u; % upper bound
end