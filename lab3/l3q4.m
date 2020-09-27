clc;
close all;
warning('off','all')

s0 = 100;
K = 100;
T = 1;
r = 0.08;
sig = 0.2;
M = [5:5:25];
compare(s0, K, T, r, sig, M);
function compare(s0, K, T, r, sig, M)
    fprintf('Using non-Markov method for American Option:\n');
    fprintf('Initial prices:\n');
    tic;
    for i=1:size(M, 2)
        [call, put] = nonmarkovOptionPricing(s0, K, T, r, sig, M(i));
        fprintf('M = %d call: %f put: %f\n', M(i), call, put);
    end
    toc;
    M = [M 50];
    fprintf('Using Markov method for American Option:\n');
    fprintf('Initial prices:\n');
    tic;
    for i=1:size(M, 2)
        [call, put] = americanOptionPricing(s0, K, T, r, sig, M(i));
        fprintf('M = %d call: %f put: %f\n', M(i), call, put);
    end
    toc;
end

function [call, put] = nonmarkovOptionPricing(s0, K, T, r, sig, M)

    [u, d] = set(T, M, sig, r);
    p  = (u-exp(r*(T/M)))/(u-d);
    %no arbitrage condition %
    if p<0||p>1
        fprintf('arbitrage detected function exiting ....%d, %d, %d, %d\n', s0, K, r, sig);
        call=0; 
        put=0;
        return;
    end
    [s, vc, vp] = finalValuesStock(s0, u, d, K, M);
    for j=1:size(vc, 2)
        vc(j) = max(vc(j), 0);
        vp(j) = max(vp(j), 0);
    end
    for j=M:-1:1
        for i = 1:2:2^j-1
            k = floor((i+1)/2);
            s(1, k) = s(1, i+1)/d;
            vc(1, k) = max((p*vc(1, i)+(1-p)*vc(1, i+1))*(exp(-r*(T/M))), s(1, k)-K);
            vp(1, k) = max((p*vp(1, i)+(1-p)*vp(1, i+1))*(exp(-r*(T/M))), K-s(1, k));
        end
    end
    call = vc(1);
    put = vp(1);
end



function [call, put] = americanOptionPricing(s0, K, T, r, sig, M)

    [u, d] = set(T, M, sig, r);
    p  = (u-exp(r*(T/M)))/(u-d);
    %no arbitrage condition %
    if p<0||p>1
        fprintf('arbitrage detected function exiting ....%d, %d, %d, %d\n', s0, K, r, sig);
        call=0; 
        put=0;
        return;
    end
    [vc, vp] = getIntrinsicValues(s0, u, d, K, M);

    for j=1:size(vc, 2)
        vc(j) = max(vc(j), 0);
        vp(j) = max(vp(j), 0);
    end
    for i=M-1:-1:0
        [gc, gp] = getIntrinsicValues(s0, u, d, K, i);
        vnc = zeros(1, i+1);
        vnp = zeros(1, i+1);
        for j=0:i
            vnc(j+1) = max((p*vc(j+1)+(1-p)*vc(j+2))*(exp(-r*(T/M))), gc(j+1));
            vnp(j+1) = max((p*vp(j+1)+(1-p)*vp(j+2))*(exp(-r*(T/M))), gp(j+1));
        end
        vc = vnc;
        vp = vnp;
    end
    
    call = vc(1);
    put = vp(1);
end

function [vc, vp] = getIntrinsicValues(s0, u, d, K, M)
    s = zeros(1, M+1);
    for j=0:M
        s(1, j+1) = s0*(u^(M-j))*(d^j);
    end 
    vc = zeros(1, M+1);
    vp = zeros(1, M+1);
    for i=0:M
        vc(1, i+1) = s(1, i+1)-K;
        vp(1, i+1) = K-s(1, i+1);
    end
end

function [s, vc, vp] = finalValuesStock(s0, u, d, K, M);
    s = zeros(1, 2^M);
    vc = zeros(1, 2^M);
    s2 = zeros(1, 2^M);
    vp = zeros(1, 2^M);
    s(1) = s0;
    for i=1:M
        for j=1:2^(i-1)
            s2(1, 2*j-1)= u*(s(1, j));
            s2(1, 2*j) = d*(s(1, j));
        end
        for j=1:2^i
            s(1, j) = s2(1, j);
        end
    end
    for i=1:2^M
        vc(1, i) = s(1, i)-K;
        vp(1, i) = K - s(1, i);
    end
end


function [u, d] = set(T, M, sig, r)
   u = exp(sig*sqrt(T/M)+(r-0.5*sig*sig)*(T/M));
   d = exp(-sig*sqrt(T/M)+(r-0.5*sig*sig)*(T/M));
end
