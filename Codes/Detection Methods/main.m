clc, clear
close all

%% Initialization
mod = 'qpsk';                                                % Modulation Name                                 
% methods_c = {'SD'};
methods_c = {'SD', 'ZF', 'MMSE', 'LRA-ZF'};

N = 3.072e5;                                                % Number of Bits 
k = 2;                                                      % Bits per Symbol
M = 2^k;                                                    % Modulation Order
Nt = 2;                                                     % Number of Transmit Antennas                                             
Nr = 2;                                                     % Number of Recieve Antennas
T = N/(k*Nt);                                               % Number of Transmission Cycles
H0 = 1;                                                     % Channel Parameter Power

% snrDB_v = 1000;
snrDB_v = 10:30;
snr_v = 10.^(snrDB_v./10);
isGray = 1;

%% Simulation
for method = methods_c
    Pe_v = [];
    for snr = snr_v
        txBit_m =  randi([0 1], N/k, k);
        [symbolIndex_v, biMatrix_m] = symbolIndexGenerator(txBit_m, N, k, isGray);
        [cons, consEnergy] = constellation(M, mod);
        z_v = cons(symbolIndex_v);
        z_m = reshape(z_v, Nt, T);
        N0 = Nt/snr;
        rxBit_m = txBit_m.*0;
        for t = 1:T
            % Transmit
            x_v = z_m(:, t);

            % Recieve
            n_v = (randn(Nr, 1) + 1j*randn(Nr, 1)) * sqrt(N0/2);
            H_m = (randn(Nr, Nt) + 1j*randn(Nr, Nt)) * sqrt(H0/2);
            y_v = H_m*x_v + n_v;

            % Process
            r_v = detector(y_v, H_m, snr, N0, Nt, Nr, cons, consEnergy, method{1}, mod);
            rxBit_m((t-1)*Nt+1:t*Nt, :) = biMatrix_m(r_v, :);
        end
        Pe_v = [Pe_v, sum(txBit_m~=rxBit_m, "all")/N];
    end
    semilogy(snrDB_v, Pe_v, 'Marker', 'x')
    hold on
end
title(['Pe_{bits}   (', mod, ', M=', num2str(M), ', Nt=', num2str(Nt), ', Nr=', num2str(Nr), ')'])
xlabel('SNR (dB)')
ylabel('Pe')
grid('on')
legend('SD', 'Zf', 'MMSE', 'LRA-ZF')

%% Repository
% helper_m = ones(M, Nr).*(0:M-1).';

%             r_v = (round(((pinv(H_m)*y_v)*sqrt(consEnergy)-(1+1j))/2)*2+(1+1j))/sqrt(consEnergy);
%             w_v = helper_m((r_v==cons).');