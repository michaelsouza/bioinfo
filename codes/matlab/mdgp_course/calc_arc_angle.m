function omega = calc_arc_angle(da,db,dc_min,dc_max)
unit  = @(v)v/norm(v);
a = [0,0,0];
b = [db 0 0]
%% calc omega
x = intersect3spheres(a,dax,b,dbx,c,dcx);
% circle center
p = a + dot(x(1)-a,unit(b-a)) * unit(b-a);
n = unit(cross(b-a,c-a));
% select x in the positive semi-space
if dot(x(1,:) - p, p + n) > 0
    x = x(1,:);
else
    x = x(2,:);
end
% circle radius
r = norm(x-p);
y = p + r * n;
q = p + r * unit(cross(y-a,b-a));
omega = acos(dot(unit(q-p),unit(x-p)));
print_angle('omega',omega);
z = rotate(omega,b-a,a,q)
print_point('z',z);