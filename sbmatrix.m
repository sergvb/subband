function m = sbmatrix( omega1, omega2, tau, N )
%SBMATRIX Summary of this function goes here
%   Detailed explanation goes here

m = zeros(N);

for k = 1:N
    for n = 1:N
        d = k - n + tau;
        m(k, n) = 1;
        if d ~= 0
            m(k, n) = -1i/d*(exp(1i*omega2*d) - exp(1i*omega1*d));
        end
    end
end

end

