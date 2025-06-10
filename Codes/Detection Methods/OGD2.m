function r_v = OGD2(y_v, H_m, cons, consEnergy, modulation)
    [yr_v] = wrapper(y_v, H_m, consEnergy, modulation);
    
    [~, m] = size(H_m);
    M = length(cons);
    switch modulation
        case 'pam'
            bounds_v = [-M/2 + 1, M/2];
        case 'qam'
            bounds_v = [-sqrt(M)/2 + 1, sqrt(M)/2];
    end

    x_v = (rand(size(H_m, 2), 1)*(bounds_v(2)-bounds_v(1)) + bounds_v(1)) +...
        1j*(rand(size(H_m, 2), 1)*(bounds_v(2)-bounds_v(1)) + bounds_v(1));
    
    cost_v = [5, 1];
    i = 0;
%     for i = 1:100*length(y_v)
    while abs(cost_v(end))>1e-4
        i = i+1;
        t = mod(i-1, size(H_m, 2))+1;
        eta = 1;
%         w_v = x_v.*H_m(t, :).';
%         w_v = w_v + 2*eta*(y_v(t)-sum(w_v))*ones(size(H_m(t, :).'));
%         x_v = w_v./H_m(t, :).';
        x_v = (x_v + eta*((y_v(t)-H_m(t, :)*x_v)/(H_m(t, :)*H_m(t, :)'))*H_m(t, :)');
%         x_v =  x_v + 2*eta*H_m.'*(y_v-H_m*x_v);

%         x_v = min(max(x_v, bounds_v(1)), bounds_v(2));
        cost_v = [cost_v, (norm(H_m*x_v-y_v, 2)^2)];
    end
    
    
    cost_v = cost_v(3:end);
    plot(cost_v)
    r_v = round(x_v);

    r_v = unwrapper(r_v, consEnergy, modulation);
    [~, r_v] = min((r_v-cons).^2, [], 2);
end

