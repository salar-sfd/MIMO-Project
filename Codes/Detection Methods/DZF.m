function [r, simTime] = DZF(y, H, cons)
    tic

    [~, Nt] = size(H);

%     func = @(x) (abs(H'*y-H'*H*x)).^2;
%     cost = zeros(Nt, Nt);

    Proj = zeros(Nt, Nt);
    yp = H(:, 1)'*y;
    h = H(:, 1)'*H;
    x = (yp/(h*h'))*h';
%     cost(1, :) = func(x);

    for i = 2:Nt
        hp_v = h' - Proj*h';
        Proj = Proj + hp_v*hp_v'/(hp_v'*hp_v);

        yp = H(:, i)'*y;
        h = H(:, i)'*H;
        dx = ((yp-h*x)/(h*h'))*h';
        dxp = dx - Proj*dx;
        dx = dxp * ((dx'*dx)/(dxp'*dxp));
        x = x + dx;     
%         cost(i, :) = func(x);
    end
    
%     for i = 1:Nt
%         subplot(Nt, 2, 2*i-1);
%         plot(cost(:, i), 'r');
%         ylabel(['c_{', num2str(i), '}', '(x)']);
%         ylim([0, max(max(cost))]);
%         xlim([1, Nt]);
%     end

%     subplot(Nt, 2, (1:Nt)*2);
%     plot(sum(cost, 2));
%     ylabel('c(x)');
%     ylim([0, max(sum(cost, 2))]);
%     xlim([1, Nt]);
%      
    [~, r] = min((x-cons).^2, [], 2);
    simTime = toc;
end





























%     y_v = [real(y_v); imag(y_v)];
%     H_m = [real(H_m), -imag(H_m); imag(H_m), real(H_m)];

%     cost_v = (norm(H_m*x_v-y_v, 2)^2);
%         i = mod(t-1, Nr)+1;
%         cost_v = [cost_v, (norm(H_m*x_v-y_v, 2)^2)];
%     plot(cost_v)

%         for j = 1:(i-1)
%             dxp_v = dxp_v - ((H_m(j, :)*dxp_v)/(H_m(j, :)*H_m(j, :)'))*H_m(j, :)';
% %             temp = (norm(H_m(i, :))*norm(H_m(j, :))/(H_m(i, :)*H_m(j, :)'));
% %             dx_v = temp*dx_v - ((y_v(i)-H_m(i, :)*x_v)/(H_m(i, :)*H_m(j, :)'))*H_m(j, :)'/temp;
%         end

%             Proj = H_m(1:(i-1), :)'*(H_m(1:(i-1), :)*H_m(1:(i-1), :)')^-1*H_m(1:(i-1), :);

% 
% function [r_v, simTime] = DZF(y_v, H_m, cons)
%     tic
%     [~, Nt] = size(H_m);
% 
%     y_v = H_m'*y_v;
%     H_m = (H_m'*H_m);
%     
%     x_v = (y_v(1)/(H_m(1, :)*H_m(1, :)'))*H_m(1, :)';
%     Proj = zeros(Nt, Nt);
% 
%     for i = 2:Nt        
%         hp_v = H_m(i-1, :)' - Proj*(H_m(i-1, :)');
%         Proj = Proj + hp_v*hp_v'/(hp_v'*hp_v);
% 
%         dx_v = ((y_v(i)-H_m(i, :)*x_v)/(H_m(i, :)*H_m(i, :)'))*H_m(i, :)';
%         dxp_v = dx_v - Proj*dx_v;
%         dx_v = dxp_v * ((dx_v'*dx_v)/(dxp_v'*dxp_v));
%         x_v = x_v + dx_v;        
% 
%     end
%     
%     [~, r_v] = min((x_v-cons).^2, [], 2);
%     simTime = toc;
% end