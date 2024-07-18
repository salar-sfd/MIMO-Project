function [H_m, U_m, u_m] = sizeReduce(H_m, U_m, u_m, k, j)
    c = round(u_m(k, j));
    H_m(:, k) = H_m(:, k)-c*H_m(:, j);
    U_m(:, k) = U_m(:, k)-c*U_m(:, j);
    for l = 1:j
        u_m(k, l) = u_m(k, l) - c*u_m(j, l);
    end
end

