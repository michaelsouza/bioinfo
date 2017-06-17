function x = intersect4spheres(a,da,b,db,c,dc,d,dd)
y = intersect3spheres(a,da,b,db,c,dc);
dev1 = abs(norm(d - y(1,:)) - dd);
dev2 = abs(norm(d - y(2,:)) - dd);
if dev1 < dev2
    x = y(1,:);
else
    x = y(2,:);
end
end