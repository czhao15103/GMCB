% Computation Effort Examination for HSGHS

%% n100p5q5 small signals, CS
% read in true b values
filename = 'n100p5q5_smallsignalsCS_ComputationComparison';
filetype = '.mat';

% extract the timing in seconds
out = [];
for j = 1:2000
    disp(j);
    istr = num2str(j);
    fileload = strcat(filename, istr, filetype);
    load(fileload)
    timing = time;
    clear time;
    converttoduration = time(timing);
    out(j) = seconds(converttoduration);
end

writematrix(transpose(out), 'hsghs_n100p5q5_smallsignalsCS_timinginsecs.csv')

%% n100p5q5 moderate sparse signals, AR1
% read in true b values
filename = 'n100p5q5_moderatesparsesignalsAR_ComputationComparison';
filetype = '.mat';

% extract the timing in seconds
out = zeros(2000, 1);
parfor (j = 1:2000,10)
    disp(j);
    istr = num2str(j);
    fileload = strcat(filename, istr, filetype);
    run = load(fileload);
    timing = getfield(run, 'time');
    converttoduration = time(timing);
    out(j) = seconds(converttoduration);
end

writematrix(out, 'hsghs_n100p5q5_moderatesparsesignalsAR_timinginsecs.csv')

%% n100p5q5 DP2002_3A_XJ2019_1
% read in true b values
filename = 'n100p5q5_DP2002_3A_XJ2019_1_ComputationComparison';
filetype = '.mat';

% extract the timing in seconds
out = zeros(2000, 1);
parfor (j = 1:2000,10)
    disp(j);
    istr = num2str(j);
    fileload = strcat(filename, istr, filetype);
    run = load(fileload);
    timing = getfield(run, 'time');
    converttoduration = time(timing);
    out(j) = seconds(converttoduration);
end

writematrix(out, 'hsghs_n100p5q5_DP2002_3A_XJ2019_1_timinginsecs.csv')