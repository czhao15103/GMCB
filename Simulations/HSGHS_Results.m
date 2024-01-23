% Result Examination for HSGHS

%% n100p5q5 small signals, CS
% read in true b values
fileloc = 'smallsignalsCS_b';
filetype = '.csv';
bfile = strcat(fileloc, filetype);
b = readmatrix(bfile);
b_cell = arrayfun(@(~)b, 1:2000, 'UniformOutput', false);

% compute true omega matrix
sigma = arrayfun(@(~)toeplitz([1 0.7 0.7 0.7 0.7]), 1:2000, 'UniformOutput', false);
omega = arrayfun(@(~)inv(sigma{1}), 1:2000, 'UniformOutput', false);

% load estimates
load('n100p5q5_smallsignalsCS.mat')

f_loss_postmean = cellfun(@(x,y,z,m)f_loss(x,y,z,m), beta_mean, b_cell, omega_mean, omega);
f_loss_steinquad = cellfun(@(x,y,z,m)f_loss(x,y,z,m), beta_quad, b_cell, omega_stein, omega);
f_loss_postmeanb = cellfun(@(x,y)norm(x - y, 'fro')^2, beta_mean, b_cell);
f_loss_postmeanomega = cellfun(@(x,y)norm(x - y, 'fro')^2, omega_mean, omega);
f_loss_quadb = cellfun(@(x,y)norm(x - y, 'fro')^2, beta_quad, b_cell);
f_loss_steinomega = cellfun(@(x,y)norm(x - y, 'fro')^2, omega_stein, omega);

loss = [f_loss_postmean; f_loss_steinquad; f_loss_postmeanb; f_loss_postmeanomega; ...
    f_loss_quadb; f_loss_steinomega]';
T = array2table(loss, 'VariableNames', {'f_loss_postmean', 'f_loss_steinquad', ...
    'f_losspostmeanb', 'f_loss_postmeanomega', 'f_loss_quadb', 'f_loss_steinomega'});
writetable(T, 'hsghs_n100p5q5_smallsignalsCS.csv')

%% n100p5q5 moderate sparse signals, AR1
% read in true b values
fileloc = 'moderatesparsesignalsAR_b';
filetype = '.csv';
bfile = strcat(fileloc, filetype);
b = readmatrix(bfile);
b_cell = arrayfun(@(~)b, 1:2000, 'UniformOutput', false);

% compute the AR(1) matrix
arsigma = zeros(5,5);
for i = 1:5
    for j = 1:5
        arsigma(i,j) = 0.7^abs(i-j);
    end
end

% compute the inverse of the AR(1) matrix as in R
arprec = zeros(5,5);
for i = 1:5
    for j = 1:5
        if i == j
            if i == 1 || i == 5
                arprec(i,j) = 1/(1 - 0.7^2);
            else
                arprec(i,j) = (1 + 0.7^2)/(1 - 0.7^2);
            end
            
        else if j == (i - 1) || j == (i + 1)
                arprec(i,j) = -0.7/(1 - 0.7^2);
            end
        end
    end
end

% compute true omega matrix
sigma = arrayfun(@(~)arsigma, 1:2000, 'UniformOutput', false);
omega = arrayfun(@(~)arprec, 1:2000, 'UniformOutput', false);

% load estimates
load('n100p5q5_moderatesparsesignalsAR.mat')

f_loss_postmean = cellfun(@(x,y,z,m)f_loss(x,y,z,m), beta_mean, b_cell, omega_mean, omega);
f_loss_steinquad = cellfun(@(x,y,z,m)f_loss(x,y,z,m), beta_quad, b_cell, omega_stein, omega);
f_loss_postmeanb = cellfun(@(x,y)norm(x - y, 'fro')^2, beta_mean, b_cell);
f_loss_postmeanomega = cellfun(@(x,y)norm(x - y, 'fro')^2, omega_mean, omega);
f_loss_quadb = cellfun(@(x,y)norm(x - y, 'fro')^2, beta_quad, b_cell);
f_loss_steinomega = cellfun(@(x,y)norm(x - y, 'fro')^2, omega_stein, omega);

loss = [f_loss_postmean; f_loss_steinquad; f_loss_postmeanb; f_loss_postmeanomega; ...
    f_loss_quadb; f_loss_steinomega]';
T = array2table(loss, 'VariableNames', {'f_loss_postmean', 'f_loss_steinquad', ...
    'f_losspostmeanb', 'f_loss_postmeanomega', 'f_loss_quadb', 'f_loss_steinomega'});
writetable(T, 'hsghs_n100p5q5_moderatesparsesignalsAR.csv')

%% n100p5q5 DP2002_3A_XJ2019_1
% read in true b values
filename = 'n100p5q5_DP2002_3A_XJ2019_1_b.csv';
b = readmatrix(filename);
b_cell = arrayfun(@(~)b, 1:2000, 'UniformOutput', false);

filename_cov = 'n100p5q5_DP2002_3A_XJ2019_1_cov.csv';
cov = readmatrix(filename_cov);
sigma = arrayfun(@(~)cov, 1:2000, 'UniformOutput', false);

filename_prec = 'n100p5q5_DP2002_3A_XJ2019_1_prec.csv';
prec = readmatrix(filename_prec);
omega = arrayfun(@(~)prec, 1:2000, 'UniformOutput', false);

% load estimates
load('n100p5q5_DP2002_3A_XJ2019_1.mat')

f_loss_postmean = cellfun(@(x,y,z,m)f_loss(x,y,z,m), beta_mean, b_cell, omega_mean, omega);
f_loss_steinquad = cellfun(@(x,y,z,m)f_loss(x,y,z,m), beta_quad, b_cell, omega_stein, omega);
f_loss_postmeanb = cellfun(@(x,y)norm(x - y, 'fro')^2, beta_mean, b_cell);
f_loss_postmeanomega = cellfun(@(x,y)norm(x - y, 'fro')^2, omega_mean, omega);
f_loss_quadb = cellfun(@(x,y)norm(x - y, 'fro')^2, beta_quad, b_cell);
f_loss_steinomega = cellfun(@(x,y)norm(x - y, 'fro')^2, omega_stein, omega);

loss = [f_loss_postmean; f_loss_steinquad; f_loss_postmeanb; f_loss_postmeanomega; ...
    f_loss_quadb; f_loss_steinomega]';
T = array2table(loss, 'VariableNames', {'f_loss_postmean', 'f_loss_steinquad', ...
    'f_losspostmeanb', 'f_loss_postmeanomega', 'f_loss_quadb', 'f_loss_steinomega'});
writetable(T, 'hsghs_n100p5q5_DP2002_3A_XJ2019_1.csv')