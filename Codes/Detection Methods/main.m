clc, clear
close all

%% Initialization
mod = 'qam';                                                % Modulation Name                                 
% methods_c = {'SD', 'ZF', 'MMSE', 'LRA-ZF'};
methods_c = {'OGD'};

N = 3.072e6;                                                % Number of Bits 
k = 2;                                                      % Bits per Symbol
M = 2^k;                                                    % Modulation Order
Nt = 4;                                                     % Number of Transmit Antennas                                             
Nr = 4;                                                     % Number of Recieve Antennas
T = N/(k*Nt);                                               % Number of Transmission Cycles
H0 = 1;                                                     % Channel Parameter Power

% snrDB_v = 20:5:40;
snrDB_v = 1000;
snr_v = 10.^(snrDB_v./10);
isGray = 1;

%% Simulation
for method = methods_c
    PeBits_v = [];
    PeSymb_v = [];
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
        PeBits_v = [PeBits_v, sum(txBit_m~=rxBit_m, "all")/N];
        PeSymb_v = [PeSymb_v, sum(sum(txBit_m~=rxBit_m, 2)~=0)/(N/k)];
    end

    subplot(2, 1, 1)
    semilogy(snrDB_v, PeBits_v, 'Marker', 'x')
    hold on
    subplot(2, 1, 2)
    semilogy(snrDB_v, PeSymb_v, 'Marker', 'x')
    hold on
end
subplot(2, 1, 1)
title(['Pe_{bits}   (', mod, ', M=', num2str(M), ', Nt=', num2str(Nt), ', Nr=', num2str(Nr), ')'])
xlabel('SNR (dB)')
ylabel('Pe')
grid('on')
legend('SD', 'ZF', 'MMSE', 'LRA-ZF')

subplot(2, 1, 2)
title(['Pe_{symb}   (', mod, ', M=', num2str(M), ', Nt=', num2str(Nt), ', Nr=', num2str(Nr), ')'])
xlabel('SNR (dB)')
ylabel('Pe')
grid('on')
legend('SD', 'ZF', 'MMSE', 'LRA-ZF')