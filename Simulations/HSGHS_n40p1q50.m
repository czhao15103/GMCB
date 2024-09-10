% read in the data
ystart = 'n40p1q50_longrange_y';
filetype = '.csv';
y = cell([1 100]);
for i = 1:100
    istr = num2str(i);
    yfile = strcat(ystart, istr, filetype);
    y{i} = readmatrix(yfile);
end

x = ones([size(y{1}, 1) 1]);

% burnin and number of iterations for sampler
burnin = 0; nmc = 5e4;
q = size(y{1}, 2);
p = size(x, 2);

% for storing the output
filenamestart = 'n40p1q50_HSGHS_rep';

% run simulation
parfor (i = 1:100,20) % maximum number of workers is 10
    % for reproducibility
    rng(i + 49)

    % run HSGHS
    start_time = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
    [beta_save,lambda_sq_save,tau_sq_save,...
    omega_save,lambda_G_sq_save,tau_G_sq_save] = HSGHS(x,y{i},burnin,nmc,eye(q));
    end_time = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
	time_taken = hours(time(between(start_time, end_time)));
    
    istr = num2str(i);
    fullfilesave = strcat(filenamestart, istr);
    parsave(fullfilesave, beta_save, omega_save, time_taken); 
end
