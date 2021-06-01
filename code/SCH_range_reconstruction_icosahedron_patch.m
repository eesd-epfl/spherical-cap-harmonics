function recursive_reconstruction = SCH_range_reconstruction_icosahedron_patch(k_range, qm_k, theta_c, eigen_table, thetas, phis, N_eps)
% Icosahedron dome reconstruction
k_min = min(k_range);
k_max = max(k_range);

C_mat = zeros(length(thetas), (k_max+1)^2);
reconstruction = zeros(length(thetas), 3);
recursive_reconstruction = zeros(length(thetas), 3, k_max + 1);
for n = k_min:k_max
    C_mat(:, n^2 + 1:n^2 + 2*n + 1) = spherical_cap_harmonic_basis(n, theta_c, ...
        eigen_table, thetas, phis, N_eps);
    for ii = n^2 + 1:n^2 + 2*n + 1
        reconstruction(:, 1) = reconstruction(:, 1) + C_mat(:, ii) .* qm_k(ii, 1);
        reconstruction(:, 2) = reconstruction(:, 2) + C_mat(:, ii) .* qm_k(ii, 2);
        reconstruction(:, 3) = reconstruction(:, 3) + C_mat(:, ii) .* qm_k(ii, 3);
    end
    recursive_reconstruction(:, :, n+1) = reconstruction;
    fprintf("\nCompleted: %3.2f %% for n = %3.0f", ((n+1)^2/(k_max + 1)^2*100), n)
end
end