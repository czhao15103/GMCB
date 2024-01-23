
function floss = f_loss(b_est, b_actual, omega_est, omega_actual)

floss = norm(b_est - b_actual, 'fro')^2 + norm(omega_est - omega_actual, 'fro')^2;

end

