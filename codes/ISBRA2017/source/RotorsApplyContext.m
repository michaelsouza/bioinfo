classdef RotorsApplyContext < handle
    properties
        J = []
        S = []
        D = []
        get_k = []
        get_j = []
    end
    methods
        function this = RotorsApplyContext(nnodes, S)
            this.J = cell(nnodes, length(S));
            this.D = cell(nnodes, nnodes, nnodes);
            this.S = S;
            this.get_j = zeros(nnodes, 1);
            this.get_k = zeros(nnodes, 1);
            for i = 1:length(S)
                this.get_k(S(i)) = i;
            end
            for i = 1:nnodes
                for j = 1:length(S)
                    if S(j) <= i
                        this.get_j(i) = S(j);
                    end
                end
            end
        end
    end
end
