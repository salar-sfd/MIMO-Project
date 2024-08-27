function r_v = SD(y_v, H_m, d, cons, consEnergy, mod)
    yr_v = wrapper(y_v, H_m, consEnergy, mod);

    y_v = [real(yr_v); imag(yr_v)];
    H_m = [real(H_m), -imag(H_m); imag(H_m), real(H_m)];

    [~, m] = size(H_m);
    dp_v = zeros(m, 1);
    z_v = zeros(m, 1);
    UBx_v = zeros(m, 1);
    x_v = zeros(m, 1);

    [Q_m, R_m] = qr(H_m);
    R_m = R_m(1:m, :);
    Q1_m = Q_m(:, 1:m);
    Q2_m = Q_m(:, m+1:end);

    yp_v = Q1_m'*y_v;

    X_m = [];
    while isempty(X_m)
        k = m;
        z_v(m) = yp_v(m);
        dp_v(m) = sqrt(d^2 - norm(Q2_m'*y_v)^2);
    
        flag = 2;
        while true
            switch flag
                case 2
                    if R_m(k, k)>=0
                        UBx_v(k) = floor((dp_v(k)+z_v(k)) / R_m(k, k));
                        x_v(k) = ceil((-dp_v(k)+z_v(k)) / R_m(k, k)) - 1;
                    else
                        UBx_v(k) = floor((-dp_v(k)+z_v(k)) / R_m(k, k));
                        x_v(k) = ceil((dp_v(k)+z_v(k)) / R_m(k, k)) - 1;
                    end
                    flag = 3;
                case 3
                    x_v(k) = x_v(k) + 1;
                    if x_v(k) <= UBx_v(k)
                        flag = 5;
                    else
                        flag = 4;
                    end
                case 4
                    k = k+1;
                    if k == m+1
                        break;
                    else
                        flag = 3;
                    end
                case 5
                    if k == 1
                        flag = 6;
                    else
                        k = k-1;
                        z_v(k) = yp_v(k) - sum(R_m(k, k+1:m).*x_v(k+1:m).');
                        dp_v(k) = sqrt(dp_v(k+1)^2 - (z_v(k+1) - R_m(k+1, k+1)*x_v(k+1))^2);
                        flag = 2;
                    end
                case 6
                    xtemp_v = unwrapper(x_v(1:m/2) + 1j*x_v(m/2+1:m), consEnergy, mod);
                    xtemp_v = [real(xtemp_v); imag(xtemp_v)];
                    if(all((min(real(cons))<=xtemp_v) & xtemp_v<=max(real(cons))))
                        X_m = [X_m, x_v];
                    end
                    flag = 3;
            end
        end
        d = d+0.2;
    end
    
    [~, indx] = min(vecnorm(y_v - H_m*X_m));
    
    r_v = X_m(1:m/2, indx) + 1j*X_m(m/2+1:m, indx);

    r_v = unwrapper(r_v, consEnergy, mod);
    [~, r_v] = min((r_v-cons).^2, [], 2);
end
