%%%% Batch file for real hi-lat data
%%%% Repeating components removed

% random seed for consistency
eta_seed = 2;
% reps
R=1;
% parameters
tau = 0.1;
B = 2;
j_min = 2;
j_max = 3;
alpha = 3;

sigma_j = B.^(-alpha/2*(j_min:j_max))*10;
% seeding random numbers for consistency
rng(eta_seed)

knots = [0 0 0 0 0.5 1 1 1 1]*pi;
%%%% where does this r come from
r = 4;

%%%% our eta?
eta = [0; randn(r, 1)];

[post_samples, Y, DA] = data_sim_AM_sigma_rep(R, eta_seed, eta, sigma_j, tau);

%%%saving the reults
filename = ['sd_data_test_1_', num2str(eta_seed), '.mat'];
save(filename)


exit