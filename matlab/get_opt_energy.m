function [u, err] = get_opt_energy(A, T, B, x0, xf, rho, S)
% Get the optimal energy given parameters

% Get u* and x*
[~, U, err] = optim_fun(A, T, B, x0, xf, rho, S);

% get energy from input
u = sum(U.^2).*(T/1000); %dt = 1/1000, T/nStep
u = norm(u);
end

