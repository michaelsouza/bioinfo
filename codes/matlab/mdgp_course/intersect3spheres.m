function x = intersect3spheres(a,dax,b,dbx,c,dcx)
n = cross(b-a,c-a);
n = n / norm(n);

A = [a-c;
     b-c;
     n];

y = [(a*a' - c*c' - dax^2 + dcx^2)/2;
     (b*b' - c*c' - dbx^2 + dcx^2)/2;
     n*a'];

p = (A \ y)';

dap = norm(p-a);
dpx = sqrt(dax^2-dap^2);
x   = real([p + dpx * n;
       p - dpx * n]);
end