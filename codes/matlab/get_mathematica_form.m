function get_mathematica_form(A)
[m,n] = size(A);
fprintf('{');
for i = 1:m
    fprintf('{');
    for j = 1:n
        fprintf('%f', A(i,j))
        if j < n
            fprintf(',')
        end
    end
    fprintf('}');
    if i < m
        fprintf(',')
    end
end
fprintf('}');
end