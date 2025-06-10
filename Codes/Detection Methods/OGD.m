function r_v = OGD(y_v, H_m, cons, consEnergy, modulation)
%     [yr_v] = wrapper(y_v, H_m, consEnergy, modulation);
    
    y_v = [real(y_v); imag(y_v)];
    H_m = [real(H_m), -imag(H_m); imag(H_m), real(H_m)];
    
    [~, m] = size(H_m);
%     M = length(cons);
%     switch modulation
%         case 'pam'
%             bounds_v = [-M/2 + 1, M/2];
%         case 'qam'
%             bounds_v = [-sqrt(M)/2 + 1, sqrt(M)/2];
%     end

%     x_v = rand(size(H_m, 2), 1)*(bounds_v(2)-bounds_v(1)) + bounds_v(1);
    x_v = zeros(size(H_m, 2), 1);
    cost_v = [5, (norm(H_m*x_v-y_v, 2)^2)];
%     i = 0;
    for i = 1:5*length(y_v)
%     while abs(cost_v(end))>1e-4
%         i = i+1;
        t = mod(i-1, size(H_m, 2))+1;
        eta = 1;

%         x_v = (x_v + eta*((y_v(t)-H_m(t, :)*x_v))*H_m(t, :).');
        x_v = (x_v + eta*((y_v(t)-H_m(t, :)*x_v)/(H_m(t, :)*H_m(t, :).'))*H_m(t, :).');
%         x_v =  x_v + 2*eta*H_m.'*(y_v-H_m*x_v);

%         x_v = min(max(x_v, bounds_v(1)), bounds_v(2));
        r_v = x_v(1:m/2) + 1j*x_v(m/2+1:m);
        [~, r_v] = min((r_v-cons).^2, [], 2);

        x_v = [real(cons(r_v)).'; imag(cons(r_v)).'];

        cost_v = [cost_v, (norm(H_m*x_v-y_v, 2)^2)];
    end
    
    cost_v = cost_v(2:end);
    plot(cost_v);
%     x_v = round(x_v);

    r_v = x_v(1:m/2) + 1j*x_v(m/2+1:m);
%     r_v = unwrapper(r_v, consEnergy, modulation);
    [~, r_v] = min((r_v-cons).^2, [], 2);
end

