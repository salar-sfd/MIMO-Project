function r_v = SD1(y_v, H_m, d, cons, consEnergy, mod)
    global k  m  R_m  x_v  z_v  yp_v  dp_v  UBx_v  X_m;
    yr_v = wrapper(y_v, H_m, consEnergy, mod);

    y_v = [real(yr_v); imag(yr_v)];
    H_m = [real(H_m), -imag(H_m); imag(H_m), real(H_m)];

    [~, m] = size(H_m);
    dp_v = zeros(m, 1);
    z_v = zeros(m, 1);
    UBx_v = zeros(m, 1);
    x_v = zeros(m, 1);

    X_m = [];

    [Q_m, R_m] = qr(H_m);
    R_m = R_m(1:m, :);
    Q1_m = Q_m(:, 1:m);
    Q2_m = Q_m(:, m+1:end);

    
    yp_v = Q1_m'*y_v;
    

    while isempty(X_m)
        k = m;
        z_v(m) = yp_v(m);
        dp_v(m) = sqrt(d^2 - norm(Q2_m'*y_v)^2);
        func2();
        d = d+0.1;
    end
    
    [~, indx] = min(vecnorm(y_v - H_m*X_m));
    
    r_v = X_m(1:m/2, indx) + 1j*X_m(m/2+1:m, indx);

    r_v = unwrapper(r_v, consEnergy, mod);
    [~, r_v] = min((r_v-cons).^2, [], 2);
end

function func2()
    global k  R_m  x_v  z_v  dp_v  UBx_v;

    if R_m(k, k)>=0
        UBx_v(k) = floor((dp_v(k)+z_v(k)) / R_m(k, k));
        x_v(k) = ceil((-dp_v(k)+z_v(k)) / R_m(k, k)) - 1;
    else
        UBx_v(k) = floor((-dp_v(k)+z_v(k)) / R_m(k, k));
        x_v(k) = ceil((dp_v(k)+z_v(k)) / R_m(k, k)) - 1;
    end
    func3();
end

function func3()
    global k  x_v  UBx_v;

    x_v(k) = x_v(k) + 1;
    if x_v(k) <= UBx_v(k)
        func5();
    else
        func4();
    end
end

function func4()
    global k  m;

    k = k+1;
    if k == m+1
        return;
    else
        func3();
    end
end

function func5()
    global k  m  R_m  x_v  z_v  yp_v  dp_v;

    if k == 1
        func6();
    else
        k = k-1;
        z_v(k) = yp_v(k) - sum(R_m(k, k+1:m).*x_v(k+1:m).');
        dp_v(k) = sqrt(dp_v(k+1)^2 - (z_v(k+1) - R_m(k+1, k+1)*x_v(k+1))^2);
        func2();
    end
end

function func6()
    global x_v  X_m;

    X_m = [X_m, x_v];
    func3();
end