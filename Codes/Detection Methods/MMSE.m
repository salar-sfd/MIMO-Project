function r_v = MMSE(y_v, H_m, snr, cons)
    Nt = size(H_m, 2);
    W_m = H_m'/(H_m*H_m' + eye(Nt)*Nt/snr);
    r_v = W_m*y_v;
    [~, r_v] = min((r_v-cons).^2, [], 2);
end