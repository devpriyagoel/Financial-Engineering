clc;
close all;
warning('off','all')

s0 = 100;
K = 100;
T = 1;
r = 0.08;
sig = 0.2;
M = [5 10 25];
MB = [5:15 20 25];
M3 = 5:5:25;
partA(s0, K, T, r, sig, M);
partB(s0, K, T, r, sig, MB);
partC(s0, K, T, r, sig, 5);
tic
l3q3(s0, K, T, r, sig, M3);
toc
function partA(s0, K, T, r, sig, M)
    fprintf('\tM\tLookback Option Price\n');
    for i=1:size(M, 2)
        loption = lookbackInitialOptionPrice(s0, K, T, r, sig, M(i));
        fprintf('%5d\t%13.5f\n', M(i), loption);
    end
end
function partB(s0, K, T, r, sig, M)
    loption = zeros(1, size(M, 2));
    for i=1:size(M, 2)
        loption(i) = lookbackInitialOptionPrice(s0, K, T, r, sig, M(i));
    end
    figure();
    plot(M, loption, 'Linewidth', 2, 'color', 'k');
    title('Comparing Values for Different M lookback option');
    xlabel('M');
    ylabel('Lookback Option Price');
end
function partC(s0, K, T, r, sig, M)
    [u, d] = set(T, M, sig, r);
    p  = (u-exp(r*(T/M)))/(u-d);
    %no arbitrage condition %
    if p<0||p>1
        fprintf('arbitrage detected function exiting ....%d, %d, %d, %d\n', s0, K, r, sig); 
        return;
    end
    v = finalValuesLookback(s0, u, d, M);
    for j=M:-1:1
        fprintf('for time t = %f\n', (T/M)*j);
        for i = 1:1:2^j
            fprintf('%8.3f ', v(i));
        end
        fprintf('\n');
        for i = 1:2:2^j-1
            k = floor((i+2)/2);
            v(1, k) = (p*v(1, i)+(1-p)*v(1, i+1))*(exp(-r*(T/M)));
        end
    end
    fprintf('for time t = 0\n');
    fprintf('%7.3f\n', v(i));
end
function l3q3(s0, K, T, r, sig, M)
    loption = zeros(1, size(M, 2));
    for i=1:size(M, 2)
        loption(i) = lookbackInitialOptionPrice(s0, K, T, r, sig, M(i));
    end
    fprintf('Initial Prices for different M non markov:\n');
    for i=1:size(M, 2)
        fprintf('M = %d : %f\n', M(i), loption(i));
    end
end
function loption = lookbackInitialOptionPrice(s0, K, T, r, sig, M)
    [u, d] = set(T, M, sig, r);
    p  = (u-exp(r*(T/M)))/(u-d);
    %no arbitrage condition %
    if p<0||p>1
        fprintf('arbitrage detected function exiting ....%d, %d, %d, %d\n', s0, K, r, sig);
        loption=0; 
        return;
    end
    v = finalValuesLookback(s0, u, d, M);
%     loption = 0.0;
%     for i=1:2^M
%         q = sum(bitget(i, 1:M+1));
%         loption = loption + (v(i))*p^(M-q)*(1-p)^q;
%     end
%     loption = loption*exp(-r*T);
    for j=M:-1:1
        for i = 1:2:2^j-1
            k = floor((i+2)/2);
            v(1, k) = (p*v(1, i)+(1-p)*v(1, i+1))*(exp(-r*(T/M)));
        end
    end
    loption = v(1);
end

function v = finalValuesLookback(s0, u, d, M)
    s = zeros(1, 2^M);
    v = zeros(1, 2^M);

    s2 = zeros(1, 2^M);
    v2 = zeros(1, 2^M);
    
    s(1) = s0;
    v(1) = s0;
    
    for i=1:M
        for j=1:2^(i-1)
            s2(1, 2*j-1)= u*(s(1, j));
            s2(1, 2*j) = d*(s(1, j));
            v2(1, 2*j-1) = max(s2(1, 2*j-1), v(1, j));
            v2(1, 2*j) = max(s2(1, 2*j), v(1, j));
        end
        for j=1:2^i
            v(1, j) = v2(1, j);
            s(1, j) = s2(1, j);
        end
    end
    for i=1:2^M
        v(i) = (v(i)-s(i));
    end
end


function [u, d] = set(T, M, sig, r)
   u = exp(sig*sqrt(T/M)+(r-0.5*sig*sig)*(T/M));
   d = exp(-sig*sqrt(T/M)+(r-0.5*sig*sig)*(T/M));
end

