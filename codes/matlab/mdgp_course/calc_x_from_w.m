function x = calc_x_from_w(d,theta,omega)
% d     : vector of distances of consecutive atoms
% theta : vector of planar angles 
% omega : vector of dihedral angles

n = length(d);
x = zeros(n,3);
ct = cos(theta(3)); st = sin(theta(3));
x(2,:) = [-d(2),0,0];
x(3,:) = [-d(2)+d(3)*ct,d(3)*st,0];
B1 = eye(4);
B2 = [-1, 0, 0, -d(2); 0, 1, 0, 0; 0, 0, -1, 0; 0, 0, 0, 1];
B3 = [-ct, -st, 0, -d(3)*ct; st, -ct, 0, d(3)*st; 0, 0, 1, 0; 0, 0, 0, 1];
B  = B1 * B2 * B3;
for k = 4:n
    ct = cos(theta(k)); st = sin(theta(k));
    cw = cos(omega(k)); sw = sin(omega(k));
    dk = d(k);
    Bk = [-ct, -st, 0, -dk*ct; st*cw, -ct*cw, -sw, dk*st*cw;st*sw, -ct*sw, cw, dk*st*sw; 0, 0, 0, 1];
    B  = B * Bk;
    v  = B * [0;0;0;1];
    x(k,:) = v(1:3);
end
end