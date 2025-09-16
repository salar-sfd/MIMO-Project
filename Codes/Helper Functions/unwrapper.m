function [r_v] = unwrapper(r_v, consEnergy, modulation)
    switch modulation
        case 'pam'
            r_v = (2*r_v - 1)/sqrt(consEnergy);
        case 'qam'
            r_v = (2*r_v - (1+1j))/sqrt(consEnergy);
        case 'qpsk'
            r_v = (r_v - (1+1j))/sqrt(consEnergy);
    end
end

