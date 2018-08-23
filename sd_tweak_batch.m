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


[post_samples, Y, theta, phi, k_theta, k_phi, Npix] = sd_tweak_sim(R, eta_seed);

%%%saving the reults
filename = ['tweak_25_sd_thin_14_', num2str(eta_seed), '.mat'];
save(filename)


exit