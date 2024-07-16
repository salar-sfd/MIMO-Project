function r_v = detector(y_v, H_m, snr, N0, Nt, Nr, cons, method)
    switch method
        case 'ZF'
            r_v = ZF(y_v, H_m, cons);
        case 'MMSE'
            r_v = MMSE(y_v, H_m, snr, Nt, cons);
    end
end