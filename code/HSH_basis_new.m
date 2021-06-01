% Optimized to work efficiently with Matlab R2020a by Mahmoud S. Shaqfa
function Y_mat = HSH_basis_new(phi, theta, degree)

% This can be faster- needs a bit of vectorisation
Y_mat = zeros(size(theta,1), (degree+1)^2);

P_cos = 2 * cos(theta) - 1; % To consider the upper part; the lower one takes 2x+1

for n = 0:degree
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
    end
end