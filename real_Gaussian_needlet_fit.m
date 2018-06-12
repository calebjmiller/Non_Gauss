function [beta_hat, f_min] = real_Gaussian_needlet_fit(negloglik1, beta_init, lb, ub, paral)

% use interior-point algorithm for large-scale problems
if ~paral
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point',...
        'Hessian', 'lbfgs', 'Display', 'iter', 'TolX', 1e-6,...
        'MaxFunEvals', 1e4);
else
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point',...
        'Hessian', 'lbfgs', 'Display', 'iter', 'TolX', 1e-6,...
        'UseParallel', 'always', 'MaxFunEvals', 1e4);
end

[beta_hat, f_min] = fmincon(negloglik1, beta_init, [], [], [], [], lb, ub,...
    [], options);

disp(['The MLE of beta is ', mat2str(round(beta_hat*1e6)/1e6)])

end