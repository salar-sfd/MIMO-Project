clc, clear
close all

%% Initialization
modulation = 'qam';                                                % Modulation Name                                 
% methods_c = {'SD', 'ZF', 'MMSE', 'LRA-ZF'};
methods_c = {'ZF', 'SZF', 'MMSE', 'LRA-ZF'};

N = 3.072e5;                                                % Number of Bits 
k = 6;                                                      % Bits per Symbol
M = 2^k;                                                    % Modulation Order
Nt = 8;                                                    % Number of Transmit Antennas                                             
Nr = 8;                                                    % Number of Recieve Antennas
T = N/(k*Nt);                                               % Number of Transmission Cycles
H0 = 1;                                                     % Channel Parameter Power

snrDB_v = 10:5:50;
% snrDB_v = [20, 40];
snr_v = 10.^(snrDB_v./10);
isGray = 1;
%% Simulation
for method = methods_c
    PeBits_v = [];
    PeSymb_v = [];
    avgSimTime = [];
    SumSimTime = 0;
    for snr = snr_v
        txBit_m =  randi([0 1], N/k, k);
        [symbolIndex_v, biMatrix_m] = symbolIndexGenerator(txBit_m, N, k, isGray);
        [cons, consEnergy] = constellation(M, modulation);
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
            [r_v, simTime] = detector(y_v, H_m, snr, N0, Nt, Nr, cons, consEnergy, method{1}, modulation);
            rxBit_m((t-1)*Nt+1:t*Nt, :) = biMatrix_m(r_v, :);
            SumSimTime = SumSimTime + simTime;
        end
        PeBits_v = [PeBits_v, sum(txBit_m~=rxBit_m, "all")/N];
        PeSymb_v = [PeSymb_v, sum(sum(txBit_m~=rxBit_m, 2)~=0)/(N/k)];
        avgSimTime = [avgSimTime, SumSimTime/T];
    end

    subplot(3, 1, 1)
    semilogy(snrDB_v, PeBits_v, 'Marker', 'x')
    hold on
    subplot(3, 1, 2)
    semilogy(snrDB_v, PeSymb_v, 'Marker', 'x')
    hold on
    subplot(3, 1, 3)
    plot(snrDB_v, avgSimTime, 'Marker', 'x')
    hold on
end
subplot(3, 1, 1)
title(['Pe_{bits}   (', modulation, ', M=', num2str(M), ', Nt=', num2str(Nt), ', Nr=', num2str(Nr), ')'])
xlabel('SNR (dB)')
ylabel('Pe')
grid('on')
legend(methods_c)

subplot(3, 1, 2)
title(['Pe_{symb}   (', modulation, ', M=', num2str(M), ', Nt=', num2str(Nt), ', Nr=', num2str(Nr), ')'])
xlabel('SNR (dB)')
ylabel('Pe')
grid('on')
legend(methods_c)

subplot(3, 1, 3)
title(['َAverage Simulation Time   (', modulation, ', M=', num2str(M), ', Nt=', num2str(Nt), ', Nr=', num2str(Nr), ')'])
xlabel('SNR (dB)')
ylabel('s')
grid('on')
legend(methods_c)