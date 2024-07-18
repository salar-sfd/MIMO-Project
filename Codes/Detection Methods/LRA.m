function r_v = LRA(y_v, H_m, cons, consEnergy, delta)
    yred_v = (y_v*sqrt(consEnergy) + (1+1j)*sum(H_m, 2))/2;
    [Hr_m, U_m] = CLLL(H_m, delta);
    W_m = pinv(Hr_m);
    z_v = round(W_m*yred_v);
%     z_v = ((((z_v)*sqrt(consEnergy)-(1+1j))/2)*2+(1+1j))/sqrt(consEnergy);
    r_v = (2*U_m*z_v - (1+1j))/sqrt(consEnergy);
    [~, r_v] = min((r_v-cons).^2, [], 2);
end

