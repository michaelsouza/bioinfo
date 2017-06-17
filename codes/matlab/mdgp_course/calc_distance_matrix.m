function D = calc_distance_matrix(x)
natoms = size(x, 1);
D = zeros(natoms);
for i = 1:natoms
    for j = (i+1):natoms
        D(i,j) = norm(x(i,:) - x(j,:));
    end
end   
end