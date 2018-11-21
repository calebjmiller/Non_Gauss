%%%% Batch file for real hi-lat data
%%%% Voroni_sphere implemented for better sampling

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


[post_samples, Y, theta, phi, k_theta, k_phi, Npix] = vs_tweak_sim(R, eta_seed);

%%%saving the results number after vs is the seed then the number of points
filename = ['vs_0_2000_tweak_25', num2str(eta_seed), '.mat'];
save(filename)


exit