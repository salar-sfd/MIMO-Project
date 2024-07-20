function r_v = SD0(y_v, H_m, d, cons, consEnergy)
    
    y_v = [real(y_v); imag(y_v)];
    H_m = [real(H_m), -imag(H_m); imag(H_m), real(H_m)];

    [m, n] = size(H_m);
    dp_v = zeros(m, 1);
    z_v = zeros(m, 1);
    UBx_v = zeros(m, 1);
    x_v = zeros(m, 1);

    X_m = [];

    [Q_m, R_m] = qr(H_m);
    R_m = R_m(1:n, :);
    Q1_m = Q_m(:, 1:n);
    Q2_m = Q_m(:, n+1:end);

    yr_v = (y_v*sqrt(consEnergy) + sum(H_m, 2))/2;

    yp_v = Q1_m'*yr_v;
    k = m;
    dp_v(m) = sqrt(d^2 - norm(Q2_m'*yr_v)^2);
    z_v(m) = yp_v(m);

    [~, ~, ~, ~, ~, X_m] = func2(k, m, R_m, x_v, z_v, yp_v, dp_v, UBx_v, X_m);
    
    [~, indx] = min(vecnorm(yr_v - H_m*X_m));
    
    r_v = X_m(1:m/2, indx) + 1j*X_m(m/2+1:m, indx);

    r_v = (2*r_v - (1+1j))/sqrt(consEnergy);
    [~, r_v] = min((r_v-cons).^2, [], 2);
end

function [k, x_v, z_v, dp_v, UBx_v, X_m] = func2(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m)
    if R_m(k, k)>=0
        UBx_v(k) = floor((dp_v(k)+z_v(k)) / R_m(k, k));
        x_v(k) = ceil((-dp_v(k)+z_v(k)) / R_m(k, k)) - 1;
    else
        UBx_v(k) = floor((-dp_v(k)+z_v(k)) / R_m(k, k));
        x_v(k) = ceil((dp_v(k)+z_v(k)) / R_m(k, k)) - 1;
    end
    [k, x_v, z_v, dp_v, UBx_v, X_m] = func3(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m);
end

function [k, x_v, z_v, dp_v, UBx_v, X_m] = func3(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m)
    x_v(k) = x_v(k) + 1;
    if x_v(k) <= UBx_v(k)
        [k, x_v, z_v, dp_v, UBx_v, X_m] = func5(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m);
    else
        [k, x_v, z_v, dp_v, UBx_v, X_m] = func4(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m);
    end
end

function [k, x_v, z_v, dp_v, UBx_v, X_m] = func4(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m)
    k = k+1;
    if k == m+1
        return;
    else
        [k, x_v, z_v, dp_v, UBx_v, X_m] = func3(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m);
    end
end

function [k, x_v, z_v, dp_v, UBx_v, X_m] = func5(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m)
    if k == 1
        [k, x_v, z_v, dp_v, UBx_v, X_m] = func6(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m);
    else
        k = k-1;
        z_v(k) = y_v(k) - sum(R_m(k, k+1:m).*x_v(k+1:m).');
        dp_v(k) = sqrt(dp_v(k+1)^2 - (z_v(k+1) - R_m(k+1, k+1)*x_v(k+1))^2);
        [k, x_v, z_v, dp_v, UBx_v, X_m] = func2(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m);
    end
end

function [k, x_v, z_v, dp_v, UBx_v, X_m] = func6(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m)
    X_m = [X_m, x_v];
    [k, x_v, z_v, dp_v, UBx_v, X_m] = func3(k, m, R_m, x_v, z_v, y_v, dp_v, UBx_v, X_m);
end