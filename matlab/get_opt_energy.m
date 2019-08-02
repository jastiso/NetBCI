function [x_final, u, err] = get_opt_energy(A, T, B, x0, xf, rho, S, tol)
% Get the optimal energy given parameters

% Get u* and x*
[X, U, err] = optim_fun(A, T, B, x0, xf, rho, S);

% get energy from input
u = sum(U.^2).*(T/1000); %dt = 1/1000, T/nStep
u = sum(u);

x_final = X(end,1:numel(x0))';
check_error(err,tol)
end

