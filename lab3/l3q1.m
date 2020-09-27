clc;
close all;
warning('off','all')

s0 = 100;
K = 100;
T = 1;
r = 0.08;
sig = 0.2;
M = 100;
s0_v = 50:2:150; %varying s0%
K_v = 50:2:150; %varying K%
r_v = 0.01:0.01:0.15; %varying r%
sig_v = 0.1:0.05:0.5; %varying sig%
M_v = 20:2:100; %varying M%
K_v2 = [95 100 105]; %for these strike prices%

[call, put] = americanOptionPricing(s0, K, T, r, sig, M);
fprintf('Initial call price: %f\n', call);
fprintf('Initial put price: %f\n', put);
sensitiveAnalysisQues1(s0, K, T, r, sig, M, s0_v, K_v, r_v, sig_v, M_v, K_v2);
function sensitiveAnalysisQues1(s0, K, T, r, sig, M, s0_v, K_v, r_v, sig_v, M_v, K_v2)

    Results(s0, K, T, r, sig, M, s0_v, 1, 'S(0)');
    Results(s0, K, T, r, sig, M, K_v, 2, 'K');
    Results(s0, K, T, r, sig, M, r_v, 3, 'r');
    Results(s0, K, T, r, sig, M, sig_v, 4, 'sigma');
    for i=1:size(K_v2, 2)
        Results(s0, K_v2(i), T, r, sig, M, M_v, 5, ['M with K = ', num2str(K_v2(i))]);
    end
end
function Results(s0, K, T, r, sig, M, v, p, s)
    call = zeros(1, size(v, 2));
    put = zeros(1, size(v, 2));
    switch p
        case 1
            for i=1:size(v, 2)
                [call(i), put(i)] = americanOptionPricing(v(i), K, T, r, sig, M);
            end
        case 2
            for i=1:size(v, 2)
                [call(i), put(i)] = americanOptionPricing(s0, v(i), T, r, sig, M);
            end
        case 3
            for i=1:size(v, 2)
                [call(i), put(i)] = americanOptionPricing(s0, K, T, v(i), sig, M);
            end
        case 4
            for i=1:size(v, 2)
                [call(i), put(i)] = americanOptionPricing(s0, K, T, r, v(i), M);
            end
        case 5
            for i=1:size(v, 2)
                [call(i), put(i)] = americanOptionPricing(s0, K, T, r, sig, v(i));
            end
            
    end
    plotting(call, put, s, v);
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

function [u, d] = set(T, M, sig, r)
   u = exp(sig*sqrt(T/M)+(r-0.5*sig*sig)*(T/M));
   d = exp(-sig*sqrt(T/M)+(r-0.5*sig*sig)*(T/M));
end

function plotting (call, put, var, v)
    figure();
    subplot(2, 1, 1);
    plot(v, call, 'Linewidth', 2, 'Color', 'm');
    title(['American Call Option Prices on varying parameter: ', var]);
    xlabel(var);
    ylabel('Call Price');
    subplot(2, 1, 2);
    plot(v, put, 'Linewidth', 2, 'Color', 'c');
    title(['American Put Option Prices on varying parameter: ', var]);
    xlabel(var);
    ylabel('Put Price');
end