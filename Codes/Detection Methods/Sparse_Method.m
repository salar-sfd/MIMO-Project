function [r, simTime] = Sparse_Method(y, H, cons)
    tic

    [~, Nt] = size(H);
    M = length(cons);
    B = zeros(Nt, Nt*M);
    for i = 1:Nt
        B(i, (i-1)*M+1:i*M) = cons.';
    end
    D = H*B;
    
    w = zeros(Nt*M, 1);
    lambda = 0.1;
    eta = 0.001;
    cost = [];
    for i = 1:5000
        grad_w = 2*D'*(D*w-y)+2*lambda*w;
        grad_lambda = w'*w;

        w = w-eta*grad_w;
        lambda = max(0, lambda-eta*grad_lambda);

        cost = [cost, (norm(y-D*w))^2+lambda*(norm(w)^2)];
    end

    plot(cost);
    r = zeros(Nt, 1);
    for i = 1:Nt
        [~, index] = max(abs(s((i-1)*M+1:i*M)));
        r(i) = cons(index);
    end
    x = B*r;
    r = x(1:Nr) + 1j*x(Nr+1:2*Nr);
    [~, r] = min((r-cons).^2, [], 2);
    simTime = toc;
end












% function [r, simTime] = SDZF(y, H, cons)
%     tic
% 
%     [~, Nt] = size(H);
% 
%     B = zeros(Nt, Nt*length(cons));
%     for i = 1:Nt
%         B(i, (i-1)*length(cons)+1:i*length(cons)) = cons.';
%     end
%     A = H*B;
%     
%     [N, ~] = size(A);
% 
% %     func = @(x) (abs(y-H*x)).^2;
% %     xx = zeros(N, N);
% %     cost = zeros(N, N);
% 
%     Proj = zeros(N, N);
%     yp = y(1);
%     a = (A(1, :)*A');
%     z = (yp/(a*a'))*a';
% 
% %     cost(1, :) = func(H'*z);
% %     xx(:, 1) = H'*z;
% 
%     for i = 2:N
%         ap_v = a' - Proj*a';
%         Proj = Proj + ap_v*ap_v'/(ap_v'*ap_v);
% 
%         yp = y(i);
%         a = (A(i, :)*A');
%         dz = ((yp-a*z)/(a*a'))*a' + z/(a*a');
%         dzp = dz - Proj*dz;
%         dz = dzp * ((dz'*dz)/(dzp'*dzp));
%         z = z + dz;     
% %         cost(i, :) = func(x);
% %         xx(:, i) = x;
%     end
%     
% %     for i = 1:Nt
% %         subplot(Nt, 2, 2*i-1);
% %         plot(cost(:, i), 'r');
% %         ylabel(['c_{', num2str(i), '}', '(x)']);
% %         ylim([0, max(max(cost))]);
% %         xlim([1, Nt]);
% %     end
% 
% %     subplot(Nt, 2, (1:Nt)*2);
% %     plot(sum(cost, 2));
% %     ylabel('c(x)');
% %     ylim([0, max(sum(cost, 2))]);
% %     xlim([1, Nt]);
% %      
% % 
% %     figure
% %     for i = 1:Nt
% %         subplot(Nt, 1, i);
% %         plot(abs(xx(i, :)), 'r');
% % %         hold on
% % %         plot(imag(xx(i, :)), 'b');
% % %         hold off
% %         ylabel(['x_{', num2str(i), '}']);
% %         xlim([1, Nt]);
% %     end
% 
%     x = B*A'*z;
%     [~, r] = min((x-cons).^2, [], 2);
%     simTime = toc;
% end