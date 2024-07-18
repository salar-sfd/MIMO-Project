function [H_m, U_m] = CLLL(H_m, delta)
    n = size(H_m, 2);
    C_v = zeros(n, 1);
    u_m = eye(n);

    for j = 1:n
        C_v(j) = H_m(:, j)'*H_m(:, j);
    end
    for j = 1:n
        for i = j+1:n
            if j==1
                u_m(i,j) = (H_m(:, j)'*H_m(:, i))/C_v(j);
            else
                u_m(i,j) = (H_m(:, j)'*H_m(:, i) - sum(conj(u_m(j, 1:j-1)).*u_m(i, 1:j-1).*C_v(1:j-1).'))/C_v(j);
            end
            C_v(i) = C_v(i) - abs(u_m(i, j))^2*C_v(j);
        end
    end


    U_m = eye(n);
    k = 2;
    while k<=n
        if (abs(real(u_m(k, k-1))) > 1/2) || (abs(imag(u_m(k, k-1))) > 1/2)
            [H_m, U_m, u_m] = sizeReduce(H_m, U_m, u_m, k, k-1);
        end

        if C_v(k) < (delta - abs(u_m(k, k-1))^2) * C_v(k-1)
            Ctemp_v = C_v;
            utemp_m = u_m;

            Ctemp_v(k-1) = C_v(k) + abs(u_m(k, k-1))^2 * C_v(k-1);
            Ctemp_v(k) = C_v(k)*(C_v(k-1)/Ctemp_v(k-1));

            utemp_m(k, k-1) = conj(u_m(k, k-1))*(C_v(k-1)/Ctemp_v(k-1));

            utemp_m(k+1:n, k-1) = u_m(k+1:n, k-1).*(utemp_m(k, k-1)) + u_m(k+1:n,k)*(C_v(k)/Ctemp_v(k-1));
            utemp_m(k+1:n, k) = u_m(k+1:n, k-1) - u_m(k+1:n, k)*u_m(k, k-1);

            utemp_m(k-1, 1:k-2) = u_m(k, 1:k-2);
            utemp_m(k, 1:k-2) = u_m(k-1, 1:k-2);

            H_m(:, [k, k-1]) = H_m(:, [k-1, k]);
            U_m(:, [k, k-1]) = U_m(:, [k-1, k]);
            C_v = Ctemp_v;
            u_m = utemp_m;
            k = max(2, k-1);
        else
            for j = k-2:-1:1
                if (abs(real(u_m(k, j))) > 1/2) || (abs(imag(u_m(k, j))) > 1/2)
                    [H_m, U_m, u_m] = sizeReduce(H_m, U_m, u_m, k, j);
                end
            end
            k = k+1;
        end
    end

end

