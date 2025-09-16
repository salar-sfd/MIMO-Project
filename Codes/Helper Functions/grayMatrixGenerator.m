function grayMatrix = grayMatrixGenerator(k)
    temp = zeros(2^k, 1);
    grayMatrix = zeros(2^k, k);
    for i = 1:2^k
        temp(i) = bitxor(i-1, bitshift(i-1, -1));
    end
    for i = 1:k
        grayMatrix(:, i) = mod(floor(temp/2^(k-i)), 2);
    end
end

