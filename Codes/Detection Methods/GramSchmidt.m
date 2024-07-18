function [Hp_m, u_m] = GramSchmidt(H_m)
    N = size(H_m, 2);
    Hp_m = H_m.*0;
    u_m = eye(N);

    Hp_m(:, 1) = H_m(:, 1);
    for i = 2:N
        u_v = conj(H_m(:, i)'*Hp_m(:, 1:i-1))./diag(Hp_m(:, 1:i-1)'*Hp_m(:, 1:i-1)).';
        u_m(1:i-1, i) = u_v;
        Hp_m(:, i) = H_m(:, i) - sum(u_v.*Hp_m(:, 1:i-1), 2);
    end
end

