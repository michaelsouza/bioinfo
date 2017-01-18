function [x,g] = rotors_apply(x,S,W)
% x: coordinates' matrices
% s: index set
% w: angles
    if ~issorted(S)
        error('The index set s must be sorted in ascending order.')
    end
    % convert to matrix
    [nrows, ncols] = size(x);
    if nrows == 3 && ncols > nrows
        x = reshape(x,3,[]);
    end
    nnodes = length(x);
    % jacobian matrix
    C = RotorsApplyContext(nnodes, S);
	
    % applying rotations
    for k = 1:length(S)
        j = S(k);
        w = W(k);
        if j < 4
            continue;
        end
        % update x
        for q = j:nnodes
            [x(:,q),C.J{q,k}] = rotate(w,x(:,q),x(:,j-1),x(:,j-2));
        end
    end
    
    % calculating jacobian
    if nargout > 1
        g = zeros(3,nnodes,length(S));
        for k = 1:length(S)
            i = S(k);
            for q = i:length(x)
                % largest rotation index that affects Xq
                j = C.get_j(q);
                g(:,q,k) = diff_x(q,j,i,C);
            end
        end
        
        % convert 3D matrix to jacobian
        % g = d x_i / d w_j
        g = reshape(g,3*nnodes,[]);
    end
end

% d Xjq / d Wi
function g = diff_x(q,j,i,C)
j_c = j;
if ~isempty(C.D{q,j_c,i})
    g = C.D{q,j,i};
    return;
end
j_q = C.get_j(q);
if j > j_q
    j = j_q;
end
if q < i || j < i || j > q
    g = [0,0,0];
else
    k = C.get_k(j);
    if i == j
        Jqk = C.J{q,k};
        u = [1;0;0;0;0;0;0;0;0;0];
        g = (Jqk * u)';
    else
        Jqk = C.J{q,k};
        jnew = max(C.S(C.S<j));
        dp = diff_x(  q,jnew,i,C);
        dy = diff_x(j-1,jnew,i,C);
        dz = diff_x(j-2,jnew,i,C);
        u = [0,dp,dy,dz];
        g = (Jqk * u')';
    end
end
C.D{q,j_c,i} = g;
end
    