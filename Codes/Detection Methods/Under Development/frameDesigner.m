function [D, u, losses] = frameDesigner(M, N, lr) 
    
%     M_list = 1:M;
    losses = [1];
%     grad = 1;
    D = normc(rand(N, M).*2-1);
    for i=1:1000
%         for m = M_list
%             R = D(:, M_list~=m) * D(:, M_list~=m).';
%             [V, L] = eig(R);
%             [~, index] = min(diag(L));
%             D(:, m) = V(:, index);
%         end
%         Gram_matrix = abs(D.' * D);
%         Gram_matrix(eye(M,'logical')) = 0;
%         loss = sum(sum(Gram_matrix))/(M*(M-1));
%         losses = [losses, loss];
%         grad = abs(losses(end-1)-losses(end));
        grad = penalized_loss_grad(D, 1);
        D = D-lr*grad;
        loss = sum(sum(abs(D.'*D - eye(M))))/(M*(M-1));
        losses = [losses, loss];
    end
    
    Gram_matrix = abs(D.' * D);
    Gram_matrix(eye(M,'logical')) = 0;
    u = max(max(Gram_matrix)); 
end

function grad = penalized_loss_grad(D, lambda)
    [~, M] = size(D);
    
    % Term A: Gram matrix part
    G = D' * D;
    grad_A = 4 * D * (G - eye(M));
    
    % Term B: penalty on norms
    norms_sq = sum(D.^2, 1);                 % 1 x M
    penalty_weights = norms_sq - 1;          % 1 x M
    grad_B = 4 * lambda * D .* penalty_weights;  % broadcasting

    % Final gradient
    grad = grad_A + grad_B;
end
