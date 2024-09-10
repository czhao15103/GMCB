% read in the data
ystart = 'n40p1q50_meancovariance_longrange_ghs_y';
filetype = '.csv';
y = cell([1 100]);
for i = 1:100
    istr = num2str(i);
    yfile = strcat(ystart, istr, filetype);
    y{i} = readmatrix(yfile);
end

% burnin and number of iterations for sampler
burnin = 0; nmc = 5e4;

% for storing the output
filenamestart = 'n40p1q50_GHS_preconly_rep';

% run simulation
parfor (i = 1:100,20) 
    % for reproducibility
    rng(i + 312)

	n = size(y{i},1);
	S = transpose(y{i}) * y{i};

    % run GHS
    start_time = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
    [omega_save,lambda_sq_save,tau_sq_save] = GHS(S,n,burnin,nmc);
    end_time = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
	time_taken = hours(time(between(start_time, end_time)));
    
    istr = num2str(i);
    fullfilesave = strcat(filenamestart, istr);
    parsave_preconly(fullfilesave, omega_save, time_taken); 
end
