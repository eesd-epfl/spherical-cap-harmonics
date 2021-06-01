function recursive_reconstruction = SCH_reconstruction_icosahedron_dome(K_max, qm_k, theta_c, eigen_table, thetas, phis, N_eps)
% function [reconstruction, recursive_reconstruction] = SCH_reconstruction_icosahedron_dome(K_max, qm_k, theta_c, eigen_table, thetas, phis)
% Icosahedron dome reconstruction

C_mat = zeros(length(thetas), (K_max+1)^2);
reconstruction = zeros(length(thetas), 3);
recursive_reconstruction = zeros(length(thetas), 3, K_max+1);
for n = 0:K_max
    C_mat(:, n^2 + 1:n^2 + 2*n + 1) = spherical_cap_harmonic_basis(n, theta_c, ...
        eigen_table, thetas, phis, N_eps);
    for ii = n^2 + 1:n^2 + 2*n + 1
        reconstruction(:, 1) = reconstruction(:, 1) + C_mat(:, ii) .* qm_k(ii, 1);
        reconstruction(:, 2) = reconstruction(:, 2) + C_mat(:, ii) .* qm_k(ii, 2);
        reconstruction(:, 3) = reconstruction(:, 3) + C_mat(:, ii) .* qm_k(ii, 3);
    end
    recursive_reconstruction(:, :, n+1) = reconstruction;
    fprintf("\nCompleted: %3.2f %% for n = %3.0f", ((n+1)^2/(K_max+1)^2*100), n)
end
end