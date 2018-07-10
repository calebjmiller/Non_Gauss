% my augmented batch file for sim_AM_sigma_rep.m
% la di da di da
% now a batch for more threads/seeds

% PATHS
%addpath(genpath('/home/cjm/Desktop/repos'))
%load('deriv_B_spline_test_v2.mat')
%%% par pool 4
%parpool(4)
%%% Keep, change n_seed = 4 was 10 R=4 was 10
eta_seed = 2;
R = 1;
n_seed = 1;

eta_est_cell = cell(1, n_seed);
sigma_j_est_cell = cell(1, n_seed);
tau_est_cell = cell(1, n_seed);

tau = 0.1;

B = 2;
j_min = 2;
j_max = 3;
alpha = 3;

sigma_j = B.^(-alpha/2*(j_min:j_max))*10;

%%%% seeding random numbers for consistency
rng(eta_seed)
knots = [0 0 0 0 0.5 1 1 1 1]*pi;
r = 4;
eta = [0; randn(r, 1)];

eta

% parfor seed = 1:n_seed
%     maxNumCompThreads(4);
%     %calling sim_AM_sigma_rep
    [post_samples, Y, DA] = irreg_theta_sim_AM_sigma_rep(R, eta_seed, eta, sigma_j, tau);
%     %storing in the cell
    %eta_est_cell{seed} = eta_est;
    %sigma_j_est_cell{seed} = sigma_j_est;
    %tau_est_cell{seed} = tau_est;
% end

%%%putting all the guesses together
%eta_est_all = cell2mat(eta_est_cell);
%sigma_j_est_all = cell2mat(sigma_j_est_cell);
%tau_est_all = cell2mat(tau_est_cell);

%%%saving the reults
filename = ['ng_irreg_theta_inner_1_', num2str(eta_seed), '.mat'];
save(filename)

%print it
%eta_est_all
%sigma_j_est_all
%tau_est_all

exit