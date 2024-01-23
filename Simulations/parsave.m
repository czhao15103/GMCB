function parsave(filename, beta_save, omega_save, time)
  save(filename, 'beta_save', 'omega_save', 'time')
end

% source: https://www.mathworks.com/matlabcentral/answers/135285-how-do-i-use-save-with-a-parfor-loop-using-parallel-computing-toolbox
