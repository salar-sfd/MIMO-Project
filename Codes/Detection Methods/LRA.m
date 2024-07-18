function r_v = LRA(y_v, H_m, snr, cons, consEnergy, delta, method)
    Nt = size(H_m, 2);
    yred_v = (y_v*sqrt(consEnergy) + (1+1j)*sum(H_m, 2))/2;
    [Hr_m, U_m] = CLLL(H_m, delta);
    switch method
        case 'ZF'
            W_m = pinv(Hr_m);
        case 'MMSE'
            W_m = Hr_m'/(Hr_m*Hr_m' + eye(Nt)*Nt/snr);
    end

    z_v = round(W_m*yred_v);
%     z_v = ((((z_v)*sqrt(consEnergy)-(1+1j))/2)*2+(1+1j))/sqrt(consEnergy);
    r_v = (2*U_m*z_v - (1+1j))/sqrt(consEnergy);
    [~, r_v] = min((r_v-cons).^2, [], 2);
end

