function x = mdgp_xinit_random(G)
    nnodes  = G.nnodes;
    max_dij = max(G.l + G.u) / 2.0;
    x = sqrt(max_dij) * rand(3, nnodes);
end