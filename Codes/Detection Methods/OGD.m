function r_v = OGD(y_v, H_m, cons, consEnergy, mod)
    [yred_v] = wrapper(y_v, H_m, consEnergy, mod);
    
    y_v = [real(yr_v); imag(yr_v)];
    H_m = [real(H_m), -imag(H_m); imag(H_m), real(H_m)];
    
end

