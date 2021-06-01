function reconstruction = SCH_reconstruction_math(K_max, qm_k, resolution, theta_c, eigen_table, N_eps)
% Mathematical reconstruction

sample_phi_range = linspace(-pi, pi, resolution);
sample_theta_range = linspace(0, theta_c-0.01, resolution);
[sample_theta, sample_phi] = meshgrid(sample_theta_range, sample_phi_range);

temp_sz = size(sample_theta); temp_sz(3) = (K_max+1)^2;
C_mat = zeros(temp_sz);
temp_sz2 = size(sample_theta); temp_sz2(3) = 3;
reconstruction = zeros(temp_sz2);
for n = 0:K_max
    C_mat(:, :, n^2 + 1:n^2 + 2*n + 1) = spherical_cap_harmonic_basis(n, theta_c, ...
        eigen_table, sample_theta, sample_phi, N_eps);
    for ii = n^2 + 1:n^2 + 2*n + 1
        reconstruction(:, :, 1) = reconstruction(:, :, 1) + C_mat(:, :, ii) .* qm_k(ii, 1);
        reconstruction(:, :, 2) = reconstruction(:, :, 2) + C_mat(:, :, ii) .* qm_k(ii, 2);
        reconstruction(:, :, 3) = reconstruction(:, :, 3) + C_mat(:, :, ii) .* qm_k(ii, 3);
    end
    fprintf("\nCompleted: %3.2f %% for n = %3.0f", ((n+1)^2/(K_max+1)^2*100), n)
end
end