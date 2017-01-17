function x = mdgp_xinit_rotors(G)
nnodes = G.nnodes;
d = sparse(G.i, G.j, (G.l + G.u) / 2, nnodes, nnodes);
% planar angles
theta = zeros(nnodes, 1);
for i = 3:nnodes
    dAB = d(i - 2, i - 1);
    dAC = d(i - 2, i);
    dBC = d(i - 1, i);
    cos_theta = (dAB^2 + dBC^2 - dAC^2)/(2*dAB*dBC);
    if abs(cos_theta) > 1
        error('The problem constraints are not compatible with xinit_rotors.')
    end
    theta(i) = acos(cos_theta);
end
% random dihedral angles
w = rand(nnodes, 1) * 2 * pi;

% set random positions
x = zeros(3, nnodes);
x(:,2) = [-d(1,2);0;0];
x(:,3) = [-d(1,2) + d(2,3) * cos(theta(3)); d(2,3) * sin(theta(3)) ; 0];
for i = 4:nnodes
    A = x(:,i-3);
    B = x(:,i-2);
    C = x(:,i-1);
    n = cross(A-B, C-B); % normal to the plane ABC
    dBC = d(i-2, i-1);
    dCD = d(i-1, i);
    D = C + dCD * (B - C) / dBC;
    % rotation on the plane ABC
    D = rotate(theta(i), D, C, C + n);
    % rotation around BC
    x(:,i) = rotate(w(i), D, B, C);
end
end
