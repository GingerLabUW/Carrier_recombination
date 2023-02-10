% Define the rate equation
function dndt = rate_equation(t, n, k_nr, k_b)
    dndt = -k_nr * n - k_b * (n^2);
end