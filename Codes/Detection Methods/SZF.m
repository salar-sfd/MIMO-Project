function [r_v, simTime] = SZF(y_v, H_m, snr, cons, consEnergy)
    tic
    [Nr, Nt] = size(H_m);
    M = length(cons);
%     samples = 2^Nr;
    samples = 256;
    N0 = Nt/snr;
    W_m = pinv(H_m);
    r_v = W_m*y_v;
    [~, r_v] = min((r_v-cons).^2, [], 2);
    r_v = cons(r_v).';

%     R_m = zeros(Nr, samples);
    N_m = (randn(Nr, samples) + 1j*randn(Nr, samples)) * sqrt(N0/2);
%     R_m = W_m*(y_v+N_m);
    R_m = (W_m*(y_v+N_m)*sqrt(consEnergy)+(1+1j))/2;
    R_m_real = max(-M/2+1, min(round(real(R_m)), M/2));
    R_m_imag = max(-M/2+1, min(round(imag(R_m)), M/2));
    R_m = R_m_real+1j*R_m_imag;
    R_m = (2*round(R_m)-(1+1j))/sqrt(consEnergy);
    
%     for i=1:samples
%         n_v = (randn(Nr, 1) + 1j*randn(Nr, 1)) * sqrt(N0);
%         R_m(:, i) = W_m*(y_v+N_m(:, i));
%         [~, R_m(:, i)] = min((R_m(:, i)-cons).^2, [], 2);
%         R_m(:, i) = cons(R_m(:, i)).';
%     end

    [uniqueVecs, ~, idx] = unique(R_m.', 'rows');
    counts = histcounts(idx, 1:max(idx)+1);
    frequentIdx = counts > floor(samples*2e-3);

    F_m = uniqueVecs(frequentIdx, :).';

    if ~ismember(r_v.', F_m.', 'rows')
        F_m = [F_m r_v]; % Append r_v as a new column
    end

    [~, idx] = min(vecnorm(y_v-H_m*F_m));

    r_v = F_m(:, idx);
    [~, r_v] = min((r_v-cons).^2, [], 2);
    simTime = toc;
end