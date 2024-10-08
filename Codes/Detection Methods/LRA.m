function r_v = LRA(y_v, H_m, snr, cons, consEnergy, delta, method, modulation)

    [yred_v] = wrapper(y_v, H_m, consEnergy, modulation);

    [Hr_m, U_m] = CLLL(H_m, delta);

    switch method
        case 'ZF'
            W_m = pinv(Hr_m);
        case 'MMSE'
            Nt = size(H_m, 2);
            W_m = Hr_m'/(Hr_m*Hr_m' + eye(Nt)*Nt/snr);
    end

    z_v = round(W_m*yred_v);
%     z_v = ((((z_v)*sqrt(consEnergy)-(1+1j))/2)*2+(1+1j))/sqrt(consEnergy);

    r_v = unwrapper(U_m*z_v, consEnergy, modulation);

    [~, r_v] = min((r_v-cons).^2, [], 2);
end

