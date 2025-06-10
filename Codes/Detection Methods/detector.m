function [r_v, simTime] = detector(y_v, H_m, snr, N0, Nt, Nr, cons, consEnergy, method, modulation)
    simTime = 0;
    switch method
        case 'ZF'
            [r_v, simTime] = ZF(y_v, H_m, cons);
        case 'MMSE'
            [r_v, simTime] = MMSE(y_v, H_m, snr, cons);
        case 'LRA-ZF'
            [r_v, simTime] = LRA(y_v, H_m, snr, cons, consEnergy, 0.75, 'ZF', modulation);
        case 'LRA-MMSE'
            [r_v, simTime] = LRA(y_v, H_m, snr, cons, consEnergy, 0.75, 'MMSE', modulation);
        case 'SD'
            [r_v, simTime] = SD(y_v, H_m, 0.7, cons, consEnergy, modulation);
        case 'OGD'
            r_v = OGD(y_v, H_m, cons, consEnergy, modulation);
        case 'DZF'
            [r_v, simTime] = DZF(y_v, H_m, cons);
        case 'DMMSE'
            [r_v, simTime] = DMMSE(y_v, H_m, snr, cons);
        case 'RHO'
            r_v = RHO(y_v, H_m, cons, consEnergy, modulation, 0.001, 4);
        case 'Sparse_Method'
            [r_v, simTime] = Sparse_Method(y_v, H_m, cons);
        case 'SZF'
            [r_v, simTime] = SZF(y_v, H_m, snr, cons, consEnergy);
        case 'SMMSE'
            [r_v, simTime] = SZF(y_v, H_m, snr, cons);
    end
end