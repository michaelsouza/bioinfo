function check_mdgp_theta()
fprintf('Checking function theta\n')
for i = 1:100
    check()
end
fprintf('   Everything is fine\n');
end
function check()
tau = rand;
xij = rand(3,1);
[~, g] = mdgp_theta(tau, xij);
f = @(y) mdgp_theta(tau, y);
g_num = gradest(f, xij)';
if norm(g-g_num) > 1E-4
    disp_vec('xij',xij);
    disp_vec('g    ', g);
    disp_vec('g_num', g_num);
    error('Failed');
end
end
