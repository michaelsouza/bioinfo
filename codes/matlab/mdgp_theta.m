function [f,g] = mdgp_theta(tau, x)
	f = sqrt(tau^2 + norm(x)^2);
	if nargout > 1
		g = x / f;
	end
end
