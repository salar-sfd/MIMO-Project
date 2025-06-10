function r_v = RHO(y_v, H_m, cons, consEnergy, modulation, sigma, d)
    [yr_v] = wrapper(y_v, H_m, consEnergy, modulation);
    
    y_v = [real(yr_v); imag(yr_v)];
    H_m = [real(H_m), -imag(H_m); imag(H_m), real(H_m)];
    [Nr, Nt] = size(H_m);
    M = length(cons);
    switch modulation
        case 'pam'
            M = length(cons);
        case 'qam'
            M = sqrt(length(cons));
    end
    
    grad = 0.*y_v;
%     for j = 1:Nr
%         for i = 1:Nt
%             for xi = -M/2+1:M/2
%             end
%         end
%     end

    w = 0;
    for n = -M/2+1:M/2
        for m = -M/2+1:M/2
            for p = -M/2+1:M/2
                for q = -M/2+1:M/2
                    w = w + exp(-(norm(y_v-H_m*[n; m; p; q])^2)/(sigma^2));
                end
            end
        end
    end
    for j = [1, 2, 3, 4]
        z = 0;
        for i = [1, 2, 3, 4]
            for n = -M/2+1:M/2
                for m = -M/2+1:M/2
                    for p = -M/2+1:M/2
                        for q = -M/2+1:M/2
                            switch i
                                case 1 
                                    t = n;
                                case 2
                                    t = m;
                                case 3 
                                    t = p;
                                case 4
                                    t = q;
                            end
                            z = z + 2*H_m(j, i)*t*exp(-(norm(y_v-H_m*[n; m; p; q])^2)/(sigma^2))/(sigma^2);
                        end
                    end
                end
            end
        end    
        grad(j) = -2*y_v(j)/(sigma^2) + z/w;
    end

    y_v = y_v + (sigma^2/2)*grad;
    
    W_m = pinv(H_m);
    r_v = round(W_m*y_v);
    r_v = r_v(1:Nt/2) + 1j*r_v(Nt/2+1:Nt);
    r_v = unwrapper(r_v, consEnergy, modulation);
    [~, r_v] = min((r_v-cons).^2, [], 2);
end

