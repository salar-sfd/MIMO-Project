function [r, simTime] = DMMSE(y, H, snr, cons)
    tic
    [~, Nt] = size(H);

    Proj = zeros(Nt, Nt);
    yp = H(:, 1)'*y;
    h = H(:, 1)'*H;
    h(1) = h(1)+Nt/snr;
    x = (yp/(h*h'))*h';

    for i = 2:Nt
        hp_v = h' - Proj*h';
        Proj = Proj + hp_v*hp_v'/(hp_v'*hp_v);

        yp = H(:, i)'*y;
        h = H(:, i)'*H;
        h(i) = h(i)+Nt/snr;
        dx = ((yp-h*x)/(h*h'))*h';
        dxp = dx - Proj*dx;
        dx = dxp * ((dx'*dx)/(dxp'*dxp));
        x = x + dx;        
    end
    
    [~, r] = min((x-cons).^2, [], 2);
    simTime = toc;
end
























%         for j = 1:(i-1)
%             dxp_v = dxp_v - ((H(j, :)*dxp_v)/(H(j, :)*H(j, :)'))*H(j, :)';
% %             temp = (norm(H(i, :))*norm(H(j, :))/(H(i, :)*H(j, :)'));
% %             dx_v = temp*dx_v - ((y(i)-H(i, :)*x_v)/(H(i, :)*H(j, :)'))*H(j, :)'/temp;
%         end

%             Proj = H(1:(i-1), :)'*(H(1:(i-1), :)*H(1:(i-1), :)')^-1*H(1:(i-1), :);