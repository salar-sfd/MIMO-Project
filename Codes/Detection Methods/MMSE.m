function [r_v, simTime] = MMSE(y_v, H_m, snr, cons)
    tic
    [Nr, Nt] = size(H_m);
    W_m = H_m'/(H_m*H_m' + eye(Nr)*Nt/snr);
    r_v = W_m*y_v;
    [~, r_v] = min((r_v-cons).^2, [], 2);
    simTime = toc;
end