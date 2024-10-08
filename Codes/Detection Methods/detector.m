function r_v = detector(y_v, H_m, snr, N0, Nt, Nr, cons, consEnergy, method, modulation)
    switch method
        case 'ZF'
            r_v = ZF(y_v, H_m, cons);
        case 'MMSE'
            r_v = MMSE(y_v, H_m, snr, cons);
        case 'LRA-ZF'
            r_v = LRA(y_v, H_m, snr, cons, consEnergy, 0.75, 'ZF', modulation);
        case 'LRA-MMSE'
            r_v = LRA(y_v, H_m, snr, cons, consEnergy, 0.75, 'MMSE', modulation);
        case 'SD'
            r_v = SD(y_v, H_m, 0.7, cons, consEnergy, modulation);
        case 'OGD'
            r_v = OGD(y_v, H_m, cons, consEnergy, modulation);
    end
end