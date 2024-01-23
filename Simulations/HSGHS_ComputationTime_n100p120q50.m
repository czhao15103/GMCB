% Computation Effort Examination for HSGHS

% n=100, p=120, q=50
% read in true b values
filename = 'n100p120q50_moderatesparsesignalsAR';
filetype = '.mat';

% extract the timing in seconds
out_sec = [];
out_hr = [];
out_day = [];
for j = 1:50
    disp(j);
    istr = num2str(j);
    fileload = strcat(filename, istr, filetype);
    load(fileload)
    timing = time;
    clear time;
    converttoduration = time(timing);
    out_sec(j) = seconds(converttoduration);
    out_hr(j) = hours(converttoduration);
    out_day(j) = day(converttoduration);
end

writematrix(transpose(out_hr), 'hsghs_n100p120q50_timinginhours.csv')
