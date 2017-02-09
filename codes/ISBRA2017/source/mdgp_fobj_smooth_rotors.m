function [f,g] = mdgp_fobj_smooth_rotors(alpha, tau, G, y, s, w)
[x,gr] = rotors_apply(y,s,w);
[f,gf] = mdgp_fobj_smooth(alpha, tau, G, x);
g = gr' * gf;
end