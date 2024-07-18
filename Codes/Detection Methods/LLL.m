function [H_m, M_m] = LLL(H_m, delta)
    Ht_m = H_m;
    N = size(H_m, 2);
    while(true)
        % Preliminaries
        [~, u_m] = GramSchmidt(H_m);
        M_m = eye(N);
    
        % Reduction
        for i = 2:N
            for j = (i-1):-1:1
                H_m(:, i) = H_m(:, i) - round(u_m(j, i)) * H_m(:, j);
                M_m(:, i) = M_m(:, i) - round(u_m(j, i)) * M_m(:, j);
            end
        end
        [Hp_m, u_m] = GramSchmidt(H_m);
    
        % Swap
        flag = 0;
        for i = 1:N-1
            if (delta*norm(Hp_m(:, i))^2 > norm(Hp_m(:, i+1) + u_m(i+1, i)*Hp_m(:, i))^2)
                H_m(:, [i, i+1]) = H_m(:, [i+1, i]);
                M_m(:, [i, i+1]) = M_m(:, [i+1, i]);
                flag = 1;
                break
            end
        end
        if flag==0
            break
        end
    end
    M_m = inv(Ht_m)*H_m;
end

