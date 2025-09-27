function [N, T] = scheduler(currentStep, maxSteps, minN, maxN, k, Nt)
    ratio = (currentStep-1)/(maxSteps-1);
    N = minN*(maxN/minN)^ratio;
    T = ceil(N/(k*Nt));
    N = T*(k*Nt);
end