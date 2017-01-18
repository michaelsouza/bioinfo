function check_rotate()
	[~,x] = mdgp_load_problem('../../instances/mdpg_1GPV_N008_EPS0.16');
    fprintf('Checking rotate function\n');
	p = x(:,4);
	y = x(:,3);
	z = x(:,2);
	theta = pi / 3;
	v = rotate(theta, p, y, z);
	disp_vec('v', v);
    for i = 1:100
        x = rand(10,1);
        [~,g] = rotate(x(1), x(2:4), x(5:7), x(8:10));
        g_num = jacobianest(@(x)rotate_handle(x), x);
        if norm(g - g_num) > 1E-4
            fprintf('|g-g_num| = %g\n', norm(g-g_num));
            fprintf('   Failed');
            keyboard;
        end
    end
    fprintf('   Everything is fine\n');
end

function x = rotate_handle(x)
theta = x(1);
p = x(2:4);
y = x(5:7);
z = x(8:10);
x = rotate(theta,p,y,z);
end