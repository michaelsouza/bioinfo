function check_mdgp_fobj_smooth_rot()
check_definition();
check_performance();
end

function check_performance()
fname = '../../instances/dmdgp_N16_D5_P0.1_S2';
[G,x] = mdgp_load_problem(fname);
tau   = rand;
alpha = rand;
ntests = 10;
profile on
for i = 1:ntests
    k = randi(G.nnodes - 3);     % number of rotations
    S = nchoosek(4:G.nnodes,k);
    j = randi(size(S, 1));
    s = S(j,:);                  % rotations' indices
    w = (2*pi) * rand(size(s));  % rotations' angles
    tic
    [~,g] = mdgp_fobj_smooth_rotors(alpha, tau, G, x, s, w);
    fprintf('(%2d/%2d) norm(g) = %8.2f (len(s) = %2d time:%3.2f)\n',i,ntests,norm(g),length(s),toc);
end
profile off
profile viewer
end

function check_definition()
fname = '../../instances/dmdgp_N8_D5_P0.1_S2';
[G,x] = mdgp_load_problem(fname);
dtol  = 0.01;
fun   = mdgp_fobj(G,x,dtol);
tau   = 0.0001;
alpha = 0.5;
s = [4,5,7];
w = (2*pi) * rand(size(s)) * 0.001;
f = mdgp_fobj_smooth_rotors(alpha, tau, G, x, s, w);
fprintf('Checking fobj_smooth\n');
fprintf('fun = %g\n',fun)
fprintf('f   = %g\n',f)
for i = 1:10
    check_definition_core(G,x);
    x = mdgp_xinit_rotors(G);
end
end

function check_definition_core(G,x)
tau   = rand;
alpha = rand;
k = randi(G.nnodes - 3);     % number of rotations
S = nchoosek(4:G.nnodes,k);
j = randi(size(S, 1));
s = S(j,:);                  % rotations' indices
w = (2*pi) * rand(size(s));  % rotations' angles

[~,g] = mdgp_fobj_smooth_rotors(alpha, tau, G, x, s, w);
f = @(w) mdgp_fobj_smooth_rotors(alpha, tau, G, x, s, w);
g_num = gradest(f,w)';
if norm(g-g_num) > 1E-4
    disp_vec('(Failed) s', s)
    fprintf('   |g-g_num| = %g\n', norm(g - g_num));
else
    disp_vec('(Passed) s', s)
end
end