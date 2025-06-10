function [r_v, simTime] = ZF(y_v, H_m, cons)
    tic
    W_m = pinv(H_m);
    r_v = W_m*y_v;
    [~, r_v] = min((r_v-cons).^2, [], 2);
    simTime = toc;
end