function f_value = real_negloglik_Gaussian_needlet(beta_all, b_mat,b_mat_deriv, samples, Npix, A_1,A_2)

n = length(samples);
beta = beta_all(1:end-1);
tau = beta_all(end);

cov_mat = real_get_cov_Gaussian_needlet(beta, b_mat,b_mat_deriv, Npix, A_1,A_2)+tau^2*eye(n);

[beta tau]

f_value = ecmnobj(reshape(samples, 1, n), zeros(n, 1), cov_mat);

end