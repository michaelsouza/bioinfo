function check_mdgp_fobj_smooth()
fname = '../../instances/dmdgp_N8_D5_P0.1_S2';
[G,x] = mdgp_load_problem(fname);
dtol  = 0.01;
fun   = mdgp_fobj(G,x,dtol);
tau   = 0.0001;
alpha = 0.5;
f = mdgp_fobj_smooth(alpha, tau, G, x);
fprintf('Checking fobj_smooth\n');
fprintf('fun = %g\n',fun)
fprintf('f   = %g\n',f)
for i = 1:10
    check(G,x);
    x = mdgp_xinit_random(G);
end
fprintf('   Everything is fine\n');
end

function check(G,xsol)
dtol  = 0.01;
fun   = mdgp_fobj(G,xsol,dtol);
tau   = rand;
alpha = rand;
[~,g] = mdgp_fobj_smooth(alpha, tau, G, xsol);
f = @(y) mdgp_fobj_smooth(alpha, tau, G, y);
x = xsol(:);
g_num = gradest(f,x)';
if norm(g-g_num) > 1E-4
    fprintf('tau   = %g\n', tau);
    fprintf('alpha = %g\n', alpha);
    fprintf('fun = %g\n',fun)
    fprintf('f   = %g\n',f)
    disp_vec('g    ', g);
    disp_vec('g_num', g_num);
    fprintf('|g-g_num| = %g\n', norm(g - g_num));
    fprintf('   Failed\n')
    keyboard
end
end