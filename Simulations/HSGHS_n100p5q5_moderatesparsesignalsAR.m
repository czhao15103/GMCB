% read in the data
ystart = 'moderatesparsesignalsAR_y';
filetype = '.csv';
y = cell([1 2000]);
for i = 1:2000
    istr = num2str(i);
    yfile = strcat(ystart, istr, filetype);
    y{i} = readmatrix(yfile);
end

xstart = 'moderatesparsesignalsAR_x.csv';
x_single = readmatrix(xstart);

% burnin and number of iterations for sampler
% based on multivariate ESS, use 2.5e4 iterations
burnin = 0; nmc = 2.5e4;
q = size(y{1}, 2);
p = size(x_single, 2);

% for storing the output
beta_mean = cell([1 2000]);
omega_mean = cell([1 2000]);
omega_stein = cell([1 2000]);
beta_quad = cell([1 2000]);

% run simulation
parfor (i = 1:2000,20) % maximum number of workers is 20
    % for reproducibility
    rng(i + 5)

    % run HSGHS
    [beta_save,lambda_sq_save,tau_sq_save,...
    omega_save,lambda_G_sq_save,tau_G_sq_save] = HSGHS(x_single,y{i},burnin,nmc,eye(q));
    
    % calculate posterior means: Bayes estimate under Frobenius loss
    beta_mean{i} = mean(beta_save, 3);
    omega_mean{i} = mean(omega_save, 3);
    
    % calculate Bayes estimate of omega under Stein's loss
    sigma_save_cell = arrayfun(@(i)inv(omega_save(:,:,i)), 1:size(omega_save,3), 'UniformOutput', false);
    sigma_save = reshape(cell2mat(sigma_save_cell), q,[],nmc); % make it into a q x q x nmc array
    sigma_mean = mean(sigma_save, 3);
    omega_stein{i} = inv(sigma_mean);
    
    % calculate Bayes estimate of B under scalar quadratic loss
    b_times_omega_cell = arrayfun(@(i)beta_save(:,:,i)*omega_save(:,:,i), 1:size(omega_save,3), 'UniformOutput', false);
    b_times_omega = reshape(cell2mat(b_times_omega_cell), p, [], nmc);
    b_times_omega_mean = mean(b_times_omega, 3);
    beta_quad{i} = b_times_omega_mean*inv(omega_mean{i});    
end

save('n100p5q5_moderatesparsesignalsAR', 'beta_mean', 'omega_mean', 'beta_quad', 'omega_stein');

