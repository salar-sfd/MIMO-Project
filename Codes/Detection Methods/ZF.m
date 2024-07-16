function r_v = ZF(y_v, H_m, cons)
    W_m = pinv(H_m);
    r_v = W_m*y_v;
    [~, r_v] = min((r_v-cons).^2, [], 2);
end