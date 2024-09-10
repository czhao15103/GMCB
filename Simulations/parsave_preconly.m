function parsave(filename, omega_save, time)
  save(filename, 'omega_save', 'time')
end

% source: https://www.mathworks.com/matlabcentral/answers/135285-how-do-i-use-save-with-a-parfor-loop-using-parallel-computing-toolbox
