function check_xinit()
fname = '../../instances/dmdgp_N8_D5_P0.1_S2';
G = mdgp_load_problem(fname);
fprintf('Generating solution using xinit_random\n')
dij_tol = 0.0;
verbose = true;
x = mdgp_xinit_random(G);
mdgp_calculate_erros(G, x, dij_tol, verbose);
fprintf('Generating solution using xinit_jjmore\n')
x = mdgp_xinit_jjmore(G);
mdgp_calculate_erros(G, x, dij_tol, verbose);
fprintf('Generating solution using xinit_rotors\n')
x = mdgp_xinit_rotors(G);
mdgp_calculate_erros(G, x, dij_tol, verbose);
end
