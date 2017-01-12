function x = xinit_random(G)
nnodes  = G.nnodes;
max_dij = max(G.l + G.u) / 2.0;
x = max_dij * nnodes * rand(nnodes, 3);
end