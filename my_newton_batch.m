% my augmented batch file for sim_AM_sigma_rep.m
% la di da di da

% PATHS
addpath(genpath('/home/cami1103/Desktop/repos/repos'))

%%% probably can delete this
%parpool(10)

%%% Keep, change n_seed = 1 was 10 R=1 was 10
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

%%% try to run for just 1 seed not on parallel so get rid of parfor
%parfor seed = 1:n_seed
    %%%remove maxNumCompThreads shouldnt be relevant
    %maxNumCompThreads(4);
    %%%calling sim_AM_sigma_rep
    % changed the call input to eta_seed
    [eta_est, sigma_j_est, tau_est] = sim_AM_sigma_rep(R, eta_seed, eta, sigma_j, tau);
    %%%storing in the cell
    % used to be entering the value at spot seed but deleted loop so made
    % it enter at 1
    eta_est_cell{1} = eta_est;
    sigma_j_est_cell{1} = sigma_j_est;
    tau_est_cell{1} = tau_est;
%end

%%%putting all the guesses together
eta_est_all = cell2mat(eta_est_cell);
sigma_j_est_all = cell2mat(sigma_j_est_cell);
tau_est_all = cell2mat(tau_est_cell);

%%%saving the reults
filename = ['sim_rep', num2str(eta_seed), '.mat'];
save(filename, 'eta_est_all', 'sigma_j_est_all', 'tau_est_all', 'eta', 'sigma_j', 'tau')

%print it
eta_est_all
sigma_j_est_all
tau_est_all

exit