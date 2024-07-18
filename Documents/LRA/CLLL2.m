function B_reduced = CLLL2(B)

% LLL algorithm of lattice reduction for complex-valued lattices

% input: B, basis matrix; output: B_reduced, reduced basis matrix

% Cong Ling, 2005

% Based on the paper published later:

% Ying Hung Gan, Cong Ling, and Wai Ho Mow, Complex lattice reduction 

% algorithm for low-complexity full-diversity MIMO detection,

% IEEE Trans. Signal Processing, vol. 57, pp. 2701-2710, July 2009. 



M =  size(B,2);         % number of columns; each column is a basis vector

delta = .75;

[Q R] = qr(B,0);          % use QR decomposition to do Gram-Schmit orthogonization

beta = abs(diag(R)).^2;      % squared length of orthogonal vectors (unnormalized)

mu = R ./ (diag(diag(R))*ones(M));                          % normalize

mu = mu.';              % to make it lower triangular

k = 2;                  % k is the stage

i_iteration = 0;

while(i_iteration < 100*M^2)                                    % loop until B cannot be reduced any more

    i_iteration = i_iteration + 1;                          % limit the maximum iteration number

    if (abs(real(mu(k,k-1))) > 0.5) | (abs(imag(mu(k,k-1))) > 0.5)   % to make it weakly reduced = abs(mu(k,k-1)) > 0.5

        [B,mu] = size_reduce_k(B,mu,k,k-1);

    end

    if(beta(k) < (delta - abs(mu(k,k-1))^2) * beta(k-1))    % swap the two if the k-1th vector is in a sense longer than the kth

        b = B(:,k);

        B(:,k) = B(:,k-1);

        B(:,k-1) = b;                                       % swap the kth column and the k-1th column

        muswap = mu(k-1,1:k-2);

        mu(k-1,1:k-2) = mu(k,1:k-2);

        mu(k,1:k-2) = muswap;                               % update mu(k-1,1:k-2) and mu(k,1:k-2)

        old_muk = mu(k+1:M,k);

        old_beta1 = beta(k-1);

        old_betak = beta(k);

        old_mu = mu(k,k-1);

        mu(k+1:M,k) = mu(k+1:M,k-1) - mu(k+1:M,k) * mu(k,k-1);

        beta(k-1) = beta(k) + abs(mu(k,k-1))^2 * beta(k-1); % update beta(k-1),beta(k), see the paper

        beta(k) = beta(k) * old_beta1 / beta(k-1);

        mu(k,k-1) = mu(k,k-1)' * old_beta1 / beta(k-1);

        mu(k+1:M,k-1) = mu(k+1:M,k-1) * mu(k,k-1) + old_muk * old_betak / beta(k-1);

        if k > 2

            k = k-1;

        end

    else

        for i = k-2 :-1: 1

            if(abs(real(mu(k,i))) > 0.5) | (abs(imag(mu(k,i))) > 0.5)

                [B,mu] = size_reduce_k(B,mu,k,i);

            end

        end

        if k < M

            k = k + 1;

        else

            B_reduced = B;

            abs(mu);

            return;

        end

    end

end

B_reduced = B;                                              % in i_iteration exceeds, simply returns B so far

'Warning: suboptimal CLLL basis'

return;



function [B,mu]=size_reduce_k(B,mu,k,j0)

% in fact on step in size reduction

eta = round(mu(k,j0));

B(:,k) = B(:,k) - eta * B(:,j0);

for i = 1 : j0-1

    mu(k,i) = mu(k,i) - eta * mu(j0,i);

end

mu(k,j0) = mu(k,j0) - eta;

return;