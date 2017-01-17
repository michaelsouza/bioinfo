function check_rotate()
	fprintf('Checking rotate function\n');
	[~,x] = mdgp_load_problem('../../instances/mdpg_1GPV_N008_EPS0.16');
	x = array2matrix(x);
	p = x(4,:);
	y = x(3,:);
	z = x(2,:);
	theta = pi / 3;
	[v, g] = rotate(theta, p, y, z);
	disp_vec('v', v);
    gt = g(:,1);
    gp = g(:,2:4);
    gy = g(:,5:7);
    gz = g(:,8:10);
	ft = @(theta) rotate(theta,p,y,z);
	disp_vec('gt    ',gt);
	disp_vec('gt_num',numdiff(ft,theta));
	fp = @(p) rotate(theta,p,y,z);
	disp_mat('gp',gp);
	disp_mat('gp_num',numdiff(fp,p));
	fy = @(y) rotate(theta,p,y,z);
	disp_mat('gy',gy);
	disp_mat('gy_num',numdiff(fy,y));
	fz = @(z) rotate(theta,p,y,z);
	disp_mat('gz',gz);
	disp_mat('gz_num',numdiff(fz,z));
end