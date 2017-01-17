function check_mdgp_phi()
fprintf('Checking function phi\n')
for i = 1:100
    check()
end
fprintf('   Everything is fine.\n')
end

function check()
alpha = rand;
tau   = rand;
y     = rand;
fun   = max(y, 0);
[~,g] = mdgp_phi(alpha, tau, y);
f = @(y) mdgp_phi(alpha, tau, y);
g_num = derivest(f,y);
if norm(g - g_num) > 1E-4
    fprintf('alpha = %g\n', alpha);
    fprintf('tau   = %g\n', tau);
    fprintf('y     = %g\n', y);
    fprintf('fun = %g\n', fun)
    fprintf('phi = %g\n', f)
    fprintf('g     = %g\n', g);
    fprintf('g_num = %g\n', g_num);
    error('Failed')
end
end