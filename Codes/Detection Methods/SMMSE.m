function [r_v, simTime] = SMMSE(y_v, H_m, snr, cons)
    tic
    [Nr, Nt] = size(H_m);
    samples = 2^Nr;
    N0 = Nt/snr;
    W_m = H_m'/(H_m*H_m' + eye(Nr)*Nt/snr);
    r_v = W_m*y_v;
    [~, r_v] = min((r_v-cons).^2, [], 2);
    r_v = cons(r_v).';

    R_m = zeros(size(H_m, 2), samples);
    for i=1:samples
        n_v = (randn(Nr, 1) + 1j*randn(Nr, 1)) * sqrt(N0/2);
        R_m(:, i) = W_m*(y_v+n_v);
        [~, R_m(:, i)] = min((R_m(:, i)-cons).^2, [], 2);
        R_m(:, i) = cons(R_m(:, i)).';
    end

    [uniqueVecs, ~, idx] = unique(R_m.', 'rows');
    counts = histcounts(idx, 1:max(idx)+1);
    frequentIdx = counts > floor(samples*0.05);

    F_m = uniqueVecs(frequentIdx, :).';

    if ~ismember(r_v.', F_m.', 'rows')
        F_m = [F_m r_v]; % Append r_v as a new column
    end

    [~, idx] = min(vecnorm(y_v-H_m*F_m));

    r_v = F_m(:, idx);
    [~, r_v] = min((r_v-cons).^2, [], 2);
    simTime = toc;
end