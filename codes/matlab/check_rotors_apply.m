function check_rotors_apply()
    fprintf('Checking rotors apply (compare with check_rotations.nb)\n')
    [~,x] = mdgp_load_problem('../../instances/mdpg_1GPV_N008_EPS0.16');
    S = [4, 5, 6];
    w = [pi/8, pi/10, pi/7];
    y = rotors_apply(x, S, w);
    disp_mat('y', y');
    
    S = 4:8;
    for k = 1:length(S)
        Sk = nchoosek(S,k);
        for j = 1:size(Sk,1)
            s = Sk(j,:);
            w = (2*pi) * rand(1,length(s));
            [~,g] = rotors_apply(x, s, w);
            f = @(w)rotors_apply_handle(x, s, w);
            g_num = jacobianest(f, w);
            if norm(g - g_num) > 1E-4
                disp_vec('(Failed) s',s);
             else
                disp_vec('(Passed) s',s);
            end
        end
    end
end

function x = rotors_apply_handle(x,s,w)
    x = rotors_apply(x,s,w);
    x = x(:);
end