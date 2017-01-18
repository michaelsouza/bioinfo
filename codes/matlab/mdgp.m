function [x,G,xsol] = mdgp()
%     fname = '../../instances/dmdgp_N8_D5_P0.1_S2';
%     fname = '../../instances/dmdgp_N16_D5_P0.1_S2';
    fname = '../../instances/dmdgp_N32_D5_P0.1_S2';
    % fname = '../../instances/dmdgp_N64_D5_P0.1_S2';
    [G, xsol] = mdgp_load_problem(fname);
    x_ini     = mdgp_xinit_rotors(G);
    % x_ini = mdgp_xinit_jjmore(G);
    
    options = sph_options();
    
    options.optim.GradObj = 'on';
    fprintf('Options\n');
    fprintf('   GradObj = %s\n', options.optim.GradObj);
    fprintf('   Plot    = %s\n\n', options.plot);
    x = sph(G,x_ini,options);
    d = rmsd(x,xsol,[],'SPH:', {'x', 'xsol'});
    fprintf('   rmsd       = %g\n\n', d);
    
    options.optim.GradObj = 'on';
    fprintf('Options\n');
    fprintf('   GradObj = %s\n', options.optim.GradObj);
    fprintf('   Plot    = %s\n', options.plot);
    x = sph_rotors(G,x_ini,options);
    d = rmsd(x,xsol,[],'SPH-ROTORS:', {'x', 'xsol'});
    fprintf('   rmsd       = %g\n\n', d); 
    
%     options.optim.GradObj = 'off';
%     fprintf('Options\n');
%     fprintf('   GradObj = %s\n', options.optim.GradObj);
%     fprintf('   Plot    = %s\n', options.plot);
%     x = sph_rotors(G,x_ini,options);
%     d = rmsd(x,xsol,[],'SPH-ROTORS:', {'x', 'xsol'});
%     fprintf('   rmsd       = %g\n\n', d); 
end