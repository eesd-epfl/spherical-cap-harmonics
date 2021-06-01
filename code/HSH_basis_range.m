function reconstruction = HSH_basis_range(phi, theta, k_range, qm_k)
k_min = min(k_range);
k_max = max(k_range);
Y_mat = zeros(size(theta,1), (k_max+1)^2);
reconstruction = zeros(length(theta), 3);

P_cos = cos(theta);

for n = k_min:k_max
    Pn = legendre(n, P_cos)';
    for m = -n:1:n
        norm = sqrt(((2*n+1) * factorial(n-abs(m)))/(2 * pi * factorial(n+abs(m))));
        if m < 0
            Y_mat(:, n^2 + n + m + 1) = Pn(:, abs(m) + 1) .* exp(1i .* abs(m) .* phi) .* norm;
        elseif m == 0
            Y_mat(:, n^2 + n + m + 1) = Pn(:, abs(m)+1) .* norm;
        elseif m > 0
            Y_mat(:, n^2 + n + m + 1) = (-1)^abs(m) .* conj(Pn(:, abs(m) + 1) .* exp(1i .* m .* phi) .* norm);
        end
        reconstruction(:, 1) = reconstruction(:, 1) + Y_mat(:, n^2 + n + m + 1) .* qm_k(n^2 + n + m + 1, 1);
        reconstruction(:, 2) = reconstruction(:, 2) + Y_mat(:, n^2 + n + m + 1) .* qm_k(n^2 + n + m + 1, 2);
        reconstruction(:, 3) = reconstruction(:, 3) + Y_mat(:, n^2 + n + m + 1) .* qm_k(n^2 + n + m + 1, 3);        
    end
end