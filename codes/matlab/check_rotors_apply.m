function check_rotors_apply()
    fprintf('Checking rotors apply (compare with check_rotations.nb)\n')
    [~,x] = mdgp_load_problem('../../instances/mdpg_1GPV_N008_EPS0.16');
    s = [4, 6, 7];
    w = [pi / 8, pi / 10, pi / 7];
    f = @(w)rotors_apply(x, s, w);
    [x,g] = rotors_apply(x, s, w);
    disp_mat('x', x);
    disp_mat('g', g);
    g_num = numdiff(f, w);
    disp_mat('g_num', g_num);
    fprintf('|g-g_num| = %g\n', norm(g(:) - g_num(:)));
end