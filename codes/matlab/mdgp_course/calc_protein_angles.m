function [omega, theta] = calc_protein_angles(x)
[n,~] = size(x);
theta = zeros(n,1);
omega = zeros(n,1);

for i = 3:n
    xa = x(i-2,:);
    xb = x(i-1,:);
    xc = x(i  ,:);
    ab = normalize(xa - xb);
    cb = normalize(xc - xb);
    theta(i) = acos(dot(ab,cb));
end

for i = 4:n
    xa = x(i-3,:);
    xb = x(i-2,:);
    xc = x(i-1,:);
    xd = x(  i,:);
    normal_abc = normalize(cross(xa - xb, xc - xb));
    normal_bcd = normalize(cross(xb - xc, xd - xc));
    omega_sign = sign(dot(cross(xd - xc,xb - xa), xc - xb));
    if omega_sign < 0
        omega(i) = 2 * pi - acos(dot(normal_abc,normal_bcd));
    else
        omega(i) = acos(dot(normal_abc,normal_bcd));
    end
end

omega = 180 * omega / pi;
theta = 180 * theta / pi;

fprintf(' id  omega  theta\n');
fprintf('%c%c%c  %c%c%c%c%c\t%c%c%c%c%c\n',175 * ones(13,1))
for i = 1:length(omega)
   fprintf('%3d  %5.1f  %5.1f\n', i, omega(i), theta(i)); 
end
end

function v = normalize(v)
    v = v / norm(v);
end