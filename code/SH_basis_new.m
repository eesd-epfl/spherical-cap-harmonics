function Y_mat = SH_basis_new(phi, theta, degree)

Y_mat = zeros(size(theta,1), (degree+1)^2);

P_cos = cos(theta); 

for n = 0:degree
    Pn = legendre(n, P_cos)';
    for m = -n:1:n
        norm = sqrt(((2*n+1) * factorial(n-abs(m)))/(4 * pi * factorial(n+abs(m))));
        if m < 0
            Y_mat(:, n^2 + n + m + 1) = Pn(:, abs(m) + 1) .* exp(1i .* abs(m) .* phi) .* norm;
        elseif m == 0
            Y_mat(:, n^2 + n + m + 1) = Pn(:, abs(m)+1) .* norm;
        elseif m > 0
            Y_mat(:, n^2 + n + m + 1) = (-1)^abs(m) .* conj(Pn(:, abs(m) + 1) .* exp(1i .* abs(m) .* phi) .* norm);
        end
    end
end