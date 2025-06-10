function r_v = RHO(y_v, H_m, cons, consEnergy, modulation, sigma, d)
    [yr_v] = wrapper(y_v, H_m, consEnergy, modulation);
    
    y_v = [real(yr_v); imag(yr_v)];
    H_m = [real(H_m), -imag(H_m); imag(H_m), real(H_m)];

    [Nr, Nt] = size(H_m);
    K = Nt/d;
    switch modulation
        case 'pam'
            M = length(cons);
        case 'qam'
            M = sqrt(length(cons));
    end
    
    indices_m = groupH(H_m, d);
    
    for i = 1:size(indices_m, 1)
        [s1_v(i), s2_m(i, :)] = calculateSums4(y_v, H_m(:, indices_m(i, :)), M, K, sigma);
    end
    grad = 0.*y_v;

%     w = prod(s1_v);
    z_v = size(Nt, 1);
    for i = 1:Nt/d
        for k = 1:d
            z_v(indices_m(i, k)) = vpa(s2_m(i, k)/s1_v(i));
        end
    end

    for j = 1:Nr
        grad(j) = -y_v(j) + sum(H_m(j, :).*z_v);
    end

    y_v = y_v + grad;
    
    W_m = pinv(H_m);
    r_v = round(W_m*y_v);
    r_v = r_v(1:Nt/2) + 1j*r_v(Nt/2+1:Nt);
    r_v = unwrapper(r_v, consEnergy, modulation);
    [~, r_v] = min((r_v-cons).^2, [], 2);
end

function indices_m = groupH(H_m, d)
    Nt = size(H_m, 2);
    indices_m = zeros(Nt/d, d);

    indices_v = 1:Nt;
    n = 1;
    while isempty(indices_v)==false
        values_v = zeros(1, length(indices_v));
        for i = 1:length(indices_v)
            values_v(i) = H_m(:, indices_v(1)).'*H_m(:, indices_v(i))/(norm(H_m(:, indices_v(1)))*norm(H_m(:, indices_v(i))));
        end
        [~, grouped_v] = maxk(abs(values_v), d);
        indices_m(n, :) = indices_v(grouped_v);
        n = n+1;
        indices_v(grouped_v) = [];
    end
end

function [s1, s2_v] = calculateSums2(y_v, H_m, M, K, sigma)
    w = 0;
    for n = -M/2+1:M/2
        for m = -M/2+1:M/2
            w = w + exp(-(norm(y_v-H_m*[n; m])^2+((1/K-1)*(y_v.'*y_v)))/(sigma^2));
        end
    end
    s1 = w;

    s2_v = zeros(1, 2);
    for i = 1:2
        z = 0;
        for n = -M/2+1:M/2
            for m = -M/2+1:M/2
                switch i
                    case 1 
                        t = n;
                    case 2
                        t = m;
                end
                z = z + t*exp(-(norm(y_v-H_m*[n; m])^2+((1/K-1)*(y_v.'*y_v)))/(sigma^2));
            end
        end
        s2_v(i) = z;
    end
end

% function [s1, s2_v] = calculateSums4(y_v, H_m, M, K, sigma)
%     w_v = zeros(1, M*M*M*M);
%     u = 1;
%     for n = -M/2+1:M/2
%         for m = -M/2+1:M/2
%             for p = -M/2+1:M/2
%                 for q = -M/2+1:M/2
%                     w_v(u) = (-(norm(y_v-H_m*[n; m; p; q])^2+((1/K-1)*(y_v.'*y_v)))/(sigma^2));
%                     u = u+1;
%                 end
%             end
%         end
%     end
%     s1 = max(w_v);
%     
%     
%     s2_v = zeros(1, 4);
%     for i = 1:4
%         z = zeros(1, M*M*M*M);
%         for n = -M/2+1:M/2
%             for m = -M/2+1:M/2
%                 for p = -M/2+1:M/2
%                     for q = -M/2+1:M/2
%                         switch i
%                             case 1 
%                                 t = n;
%                             case 2
%                                 t = m;
%                             case 3
%                                 t = p;
%                             case 4
%                                 t = q;
%                         end
%                         z_v(u) = t*exp(-(norm(y_v-H_m*[n; m; p; q])^2+((1/K-1)*(y_v.'*y_v)))/(sigma^2));
%                     end
%                 end
%             end
%         end
%         s2_v(i) = max(z_v);
%     end
% end

% function [s1, s2_v] = calculateSums4(y_v, H_m, M, K, sigma)
%     w = 0;
%     for n = -M/2+1:M/2
%         for m = -M/2+1:M/2
%             for p = -M/2+1:M/2
%                 for q = -M/2+1:M/2
%                     w = w + exp(-(norm(y_v-H_m*[n; m; p; q])^2+((1/K-1)*(y_v.'*y_v))/1.5)/(sigma^2));
%                 end
%             end
%         end
%     end
%     s1 = w;
% 
%     s2_v = zeros(1, 4);
%     for i = 1:4
%         z = 0;
%         for n = -M/2+1:M/2
%             for m = -M/2+1:M/2
%                 for p = -M/2+1:M/2
%                     for q = -M/2+1:M/2
%                         switch i
%                             case 1 
%                                 t = n;
%                             case 2
%                                 t = m;
%                             case 3
%                                 t = p;
%                             case 4
%                                 t = q;
%                         end
%                         z = z + t*exp(-(norm(y_v-H_m*[n; m; p; q])^2+((1/K-1)*(y_v.'*y_v))/1.5)/(sigma^2));
%                     end
%                 end
%             end
%         end
%         s2_v(i) = z;
%     end
% end

function [s1, s2_v] = calculateSums4(y_v, H_m, M, K, sigma)
    syms n m p q real;
    M = sym(M);         
    K = sym(K);         
    sigma = sym(sigma); 

    w = sym(0);
    
    for n = -M/2+1:M/2
        for m = -M/2+1:M/2
            for p = -M/2+1:M/2
                for q = -M/2+1:M/2
                    expr = exp(-(norm(y_v - H_m*[n; m; p; q])^2 + ((1/K-1)*(y_v.'*y_v))/1.5) / (sigma^2));
                    w = w + expr;  
                end
            end
        end
    end
    
    s1 = simplify(w);  

    s2_v = sym(zeros(1, 4));
    
    for i = 1:4
        z = sym(0);
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
                        
                        expr = t * exp(-(norm(y_v - H_m*[n; m; p; q])^2 + ((1/K-1)*(y_v.'*y_v))/1.5) / (sigma^2));
                        z = z + expr;
                    end
                end
            end
        end
        s2_v(i) = simplify(z); 
    end
end
