function [rmsd,x,y] = calc_rmsd(x,y)
%Root Mean Square Deviation
xm = mean(x);
ym = mean(y);
natoms = size(x,1);
for i = 1:natoms
    x(i,:) = x(i,:) - xm;
    y(i,:) = y(i,:) - ym;
end
[u,~,v] = svd(x'*y);
q = u*v';
x = x * q;
rmsd = norm(x-y,'fro') / sqrt(natoms);
end