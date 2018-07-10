function [post_samples, Y, DA] = irreg_theta_sim_AM_sigma_rep(R, seed, eta, sigma_j, tau)
% simulation study
% generate non-stationary std function using random seed eta_seed

nu = 4;

% the grid
%addpath(genpath('/home/cjm/Desktop/repos'))
%addpath(genpath('/home/cami1103/Desktop/repos/repos'))
B = 2;
Nside = 2^3;
tp = pix2ang(Nside, 'nest', false);
N = length(tp);
theta0 = zeros(N, 1);
phi0 = zeros(N, 1);
for i = 1:N
    theta0(i) = tp{i}(1);
    phi0(i) = tp{i}(2);
end

j_min = 2;
j_max = 3;
% make sure these match in the R file for derivatives
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
r = 4;

eta_est = zeros(r+1, R);
sigma_j_est = zeros(j_max-j_min+1, R);
tau_est = zeros(1, R);

rng(seed)
for rep = 1:R
    rep
    % perturbation possibly same grid
    theta = theta0+randn(N, 1)*pi/10;
    theta(theta<0) = theta(theta<0)+pi;
    theta(theta>pi) = theta(theta>pi)-pi;
    phi = phi0+randn(N, 1)*2*pi/10;
    phi(phi<0) = phi(phi<0)+2*pi;
    phi(phi>2*pi) = phi(phi>2*pi)-2*pi;
    %%%% lets remove a sector
   
  
    ind = theta>pi/2;
    phi = phi(ind);
    theta=theta(ind);
    
    %max(theta)
    %min(theta)
    %mean(theta)
    %save('irreg_theta_inner.mat','theta')
    %pause
    % load theta phi k_theta k_phi
  
    % design matrix A
    k_theta = pi/4;
    k_phi = pi/4;
    [Npix, ~, A_1,A_2] = get_A_full(B, j_min, j_max, theta, phi, k_theta, k_phi);
    M = size(A_1, 2); 
    N = size(A_1,1);
    % non-stationary variance function
    [b_mat, ~] = bspline_basismatrix(4, knots, theta);
    b_mat(:, 1) = 1;
    % g
    std_vec = exp(b_mat*eta);
    %g'
    load('irreg_theta_inner_deriv_B_spline.mat')
    b_mat_deriv = bS;
    % first column of b_mat_deriv to zeros consistent with the ones in b_mat
    b_mat_deriv(:,1)=0;
    std_vec_deriv = b_mat_deriv*eta.*std_vec;
    % form some fakish data notice for loop and mult rather than diagonal
    % matrix multiplication to save time
    DA = zeros(N, M);
    for i = 1:N
        DA(i, :) = std_vec(i)*A_1(i, :) - std_vec_deriv(i)*A_2(i,:);
    end

    c = zeros(M, 1);
    st = 1;
    for j = j_min:j_max
        index_j = j-j_min+1;
        range = st:st+Npix(index_j)-1;
        c(range) = sigma_j(index_j)*trnd(nu, Npix(index_j), 1);
        st = st+Npix(index_j);
    end
    
    %%%%% Change this to load fake_data
    Y = DA*c+randn(N, 1)*tau;
    
    % get init values for MCMC
    beta_init = [zeros(1, r) 1 0.1^2 1e-2];
    negloglik1 = @(beta_all) real_negloglik_Gaussian_needlet(beta_all, b_mat,b_mat_deriv, Y, Npix, A_1,A_2);

    lb = [-10*ones(1, r) 0 0 1e-3];
    ub = [10*ones(1, r) 20 20 Inf];

    [beta_hat, f_min] = real_Gaussian_needlet_fit(negloglik1, beta_init, lb, ub, false);

    % init
    % c
    c_init = zeros(M, 1);
    % V
    V_inv_init = ones(M, 1); 
    % sigma_j_sq
    sigma_j_sq_init = beta_hat(end-2:end-1)'./(nu/(nu-2));
    % eta
    eta_init = [0; beta_hat(1:r)'];
    % pri_sig of eta
    tau_eta_sq = 1e2;
    % tau
    tau_init = beta_hat(end);
    tau_sq_inv_init = 1/tau_init^2;
    % tuning parameters
    mu_init = zeros(r, 1);
    Sigma_init = eye(r);
    lambda = 0.05;
    % the number of MCMC iterations
    T = 4e5;
    % the length of the burn-in period
    burn_in = 2e5;
    % the length of the thinning interval
    thin = 200;
    % the length of the interval to report progress
    n_report = 1e3;

    model = struct('A_1', A_1, 'A_2', A_2, 'b_mat', b_mat, 'b_mat_deriv', b_mat_deriv, 'nu', nu);

    data = struct('Y', Y, 'Npix', Npix);

    params = struct('c', c_init, 'V', V_inv_init, 'sigma_j_sq', sigma_j_sq_init,...
        'eta', eta_init, 'tau_eta_sq', tau_eta_sq, 'tau', tau_sq_inv_init);

    tuning = struct('mu', mu_init, 'Sigma', Sigma_init, 'lambda', lambda);

    options = struct('T', T, 'burn_in', burn_in, 'thin', thin, 'n_report', n_report, 'save', false);

    post_samples = real_Gibbs_sampler_AM_reparam(model, data, params, tuning, options);
    
    %eta_est(:, rep) = mean(post_samples.eta, 2);
    %sigma_j_est(:, rep) = mean(sqrt(post_samples.sigma_j_sq), 2);
    %tau_est(rep) = mean(1./sqrt(post_samples.tau_sq_inv));
    
end

end