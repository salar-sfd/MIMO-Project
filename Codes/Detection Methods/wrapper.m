function [yred_v] = wrapper(y_v, H_m, consEnergy, mod)
    switch mod
        case 'pam'
            yred_v = (y_v*sqrt(consEnergy) + sum(H_m, 2))/2;
        case 'qam'
            yred_v = (y_v*sqrt(consEnergy) + (1+1j)*sum(H_m, 2))/2;
        case 'qpsk'
            yred_v = (y_v*sqrt(consEnergy) + (1+1j)*sum(H_m, 2));
    end
end

