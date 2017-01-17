function [x,g] = rotors_apply(x,S,W)
% x: coordinates' matrices
% s: index set
% w: angles
    if ~issorted(S)
        error('The index set s must be sorted in ascending order.')
    end
    [nrows, ncols] = size(x);
    if ncols > 1
        ismatrix = false;
        x = array2matrix(x);
        [nrows, ncols] = size(x);
    end
    % jacobian matrix
    J = zeros(nrows,3,10,length(S));
	
    % applying rotations
    for k = 1:length(S)
        j = S(k);
        w = W(k);
        if j < 4
            continue;
        end
        % update x
        for q = j:length(x)
            [x(q,:),J(q,:,:,k)] = rotate(w,x(q,:),x(j-1,:),x(j-2,:));
        end
    end
    
    g = zeros(nrows,ncols,length(S));
    for k = 1:length(S)
        i = S(k);
        for q = i:length(x)
            % largest rotation index that affects Xq
            j = max(S(S<=q));
            g(q,:,k) = diff_x(q,j,i,S,J);
        end
    end
    
    % convert 3D matrix to jacobian
    % g = d x_i / d w_j
    if ~ismatrix
        g_vec = zeros(size(g(:)));
        g = g_vec;
    end
end

function g = diff_x(q,j,i,S,J)
% d Xjq / d Wi
    if q < i || j < i
        g = [0,0,0];
    elseif i == j
        [~,k] = ismember(i, S); 
        Jqk = reshape(J(q,:,:,k), [3, 10]);
        u = [1;0;0;0;0;0;0;0;0;0];
        g = (Jqk * u)';
    else
        jnew = max(S(S<j));
        dp = diff_x(  q,jnew,i,S,J);
        dy = diff_x(j-1,jnew,i,S,J);
        dz = diff_x(j-2,jnew,i,S,J);
        [~,k] = ismember(j,S);
        Jqk = reshape(J(q,:,:,k),[3,10]);
        u = [0,dp,dy,dz];
        g = (Jqk * u')';
    end
end