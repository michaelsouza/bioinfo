function [f,g] = mdgp_phi(alpha, tau, x)
	f = alpha * x + sqrt((alpha * x).^2 + tau.^2);
	if nargout > 1
		g = alpha * f ./ (f - alpha * x);
	end
end
