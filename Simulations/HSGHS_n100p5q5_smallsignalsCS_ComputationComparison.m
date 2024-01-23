% read in the data
ystart = 'smallsignalsCS_y';
filetype = '.csv';
y = cell([1 2000]);
for i = 1:2000
    istr = num2str(i);
    yfile = strcat(ystart, istr, filetype);
    y{i} = readmatrix(yfile);
end

xstart = 'smallsignalsCS_x.csv';
x_single = readmatrix(xstart);

% burnin and number of iterations for sampler
burnin = 0; nmc = 1e5;
q = size(y{1}, 2);
p = size(x_single, 2);

% for storing the output
filenamestart = 'n100p5q5_smallsignalsCS_ComputationComparison';

% run simulation
parfor (i = 1:2000,20) % maximum number of workers is 20
    % for reproducibility
    rng(i + 5)

    % run HSGHS
	start_time = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
    [beta_save,lambda_sq_save,tau_sq_save,...
    omega_save,lambda_G_sq_save,tau_G_sq_save] = HSGHS(x_single,y{i},burnin,nmc,eye(q));
	end_time = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
	time_taken = between(start_time, end_time);
    
    % covariance matrix can be computed later    
    istr = num2str(i);
    fullfilesave = strcat(filenamestart, istr);
    parsave(fullfilesave, beta_save, omega_save, time_taken);
end
