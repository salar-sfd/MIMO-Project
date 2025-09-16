function [symbolIndex_v, grayMatrix] = symbolIndexGenerator(txBit_m, N, k, isGray)
    if isGray
        grayMatrix = grayMatrixGenerator(k);
        symbolIndex_v = zeros(N/k, 1);
        for i = 1:(N/k)
            symbolIndex_v(i) = find(all(grayMatrix == txBit_m(i, :), 2));
        end
    else
        symbolIndex_v = bi2de(txBit_m, 'left-msb') + 1;
    end
end