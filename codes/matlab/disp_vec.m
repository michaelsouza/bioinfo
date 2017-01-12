function disp_vec(name, x)
    fprintf('%s =', name)
    for i = 1:length(x)
        fprintf(' % g', x(i))
    end
    fprintf('\n');
end