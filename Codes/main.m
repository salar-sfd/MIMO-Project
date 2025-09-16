clc, clear
close all
t_id = tic;
%% Initialization
modulation = 'qam';                                         % Modulation Name                                 
methods_c = {'SD', 'ZF', 'MMSE', 'LRA-ZF'};

N = 3.072e3;                                                % Number of Bits 
k = 6;                                                      % Bits per Symbol
M = 2^k;                                                    % Modulation Order
Nt = 64;                                                    % Number of Transmit Antennas                                             
Nr = 64;                                                    % Number of Recieve Antennas
T = N/(k*Nt);                                               % Number of Transmission Cycles
H0 = 1;                                                     % Channel Parameter Power

isGray = 1;
snrDB_v = 10:10:40;
snr_v = 10.^(snrDB_v./10);
[cons, consEnergy] = constellation(M, modulation);

%% Simulation
for method = methods_c
    PeBits_v = [];
    PeSymb_v = [];
    avgSimTime = [];
    detector = makeDetector(method{1}, cons, consEnergy, modulation);
    for snr = snr_v
        SumSimTime = 0;
        txBit_m =  randi([0 1], N/k, k);
        [txSymbolIndex_v, biMatrix_m] = symbolIndexGenerator(txBit_m, N, k, isGray);
        z_v = cons(txSymbolIndex_v);
        z_m = reshape(z_v, Nt, T);
        N0 = Nt/snr;
        rxSymbolIndex_m = zeros(Nt, T);
        parfor t = 1:T
            % Transmit
            x_v = z_m(:, t);

            % Recieve
            n_v = (randn(Nr, 1) + 1j*randn(Nr, 1)) * sqrt(N0/2);
            H_m = (randn(Nr, Nt) + 1j*randn(Nr, Nt)) * sqrt(H0/2);
            y_v = H_m*x_v + n_v;

            % Process
            [rxSymbolIndex_m(:, t), simTime] = detector(y_v, H_m, snr, N0, Nt, Nr);
            SumSimTime = SumSimTime + simTime;
        end
        for t = 1:T
            rxBit_m((t-1)*Nt+1:t*Nt, :) = biMatrix_m(rxSymbolIndex_m(:, t), :);
        end
        PeBits_v = [PeBits_v, sum(txBit_m~=rxBit_m, "all")/N];
        PeSymb_v = [PeSymb_v, sum(sum(txBit_m~=rxBit_m, 2)~=0)/(N/k)];
        avgSimTime = [avgSimTime, SumSimTime/T];
    end

    subplot(4, 2, [1, 3, 5])
    semilogy(snrDB_v, PeBits_v, 'Marker', 'x')
    hold on
    subplot(4, 2, [2, 4, 6])
    semilogy(snrDB_v, PeSymb_v, 'Marker', 'x')
    hold on
    subplot(4, 1, 4)
    plot(snrDB_v, avgSimTime, 'Marker', 'x')
    hold on
end
subplot(4, 2, [1, 3, 5])
title(['Pe_{bits}   (', modulation, ', M=', num2str(M), ', Nt=', num2str(Nt), ', Nr=', num2str(Nr), ')'])
xlabel('SNR (dB)')
ylabel('Pe')
grid('on')
legend(methods_c, 'Location', 'southwest')

subplot(4, 2, [2, 4, 6])
title(['Pe_{symb}   (', modulation, ', M=', num2str(M), ', Nt=', num2str(Nt), ', Nr=', num2str(Nr), ')'])
xlabel('SNR (dB)')
ylabel('Pe')
grid('on')
legend(methods_c, 'Location', 'southwest')

subplot(4, 1, 4)
title(['َAverage Simulation Time   (', modulation, ', M=', num2str(M), ', Nt=', num2str(Nt), ', Nr=', num2str(Nr), ')'])
xlabel('SNR (dB)')
ylabel('s')
grid('on')
legend(methods_c)
toc(t_id)