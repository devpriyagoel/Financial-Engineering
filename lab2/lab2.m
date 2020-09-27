clc;
close all;
warning('off','all')

s0 = 100;
K = 100;
T = 1;
r = 0.08;
sig = 0.2;
M = 100;
M2=10;
s0_v = 50:2:150; %varying s0%
K_v = 50:2:150; %varying K%
r_v = 0.01:0.01:0.15; %varying r%
sig_v = 0.1:0.05:0.5; %varying sig%
M_v = 20:1:100; %varying M%
M_v2 = 1:1:20;
K_v2 = [95 100 105];

sensitiveAnalysis(s0, K, T, r, sig, M, M2, s0_v, K_v, r_v, sig_v, M_v, M_v2, K_v2);
function sensitiveAnalysis(s0, K, T, r, sig, M, M2, s0_v, K_v, r_v, sig_v, M_v, M_v2, K_v2) 
    for f=1:2
        Results(s0, K, T, r, sig, M, M2, f, s0_v, 1, 'S(0)');
        Results(s0, K, T, r, sig, M, M2, f, K_v, 2, 'K');
        Results(s0, K, T, r, sig, M, M2, f, r_v, 3, 'r');
        Results(s0, K, T, r, sig, M, M2, f, sig_v, 4, 'sigma');
    end
    varyM(s0, T, r, sig, M_v, M_v2, K_v2);
end
function Results(s0, K, T, r, sig, M, M2, f, v, p, s)
    call = zeros(1, size(v, 2));
    put = zeros(1, size(v, 2));
    asian = zeros(1, size(v, 2));
    switch p
        case 1
            for i=1:size(v, 2)
                [call(i), put(i)] = optionpricing(v(i), K, T, r, sig, M, f);
                asian(i) = AsianOptionPrices(v(i), K, T, r, sig, M2, f);
            end
        case 2
            for i=1:size(v, 2)
                [call(i), put(i)] = optionpricing(s0, v(i), T, r, sig, M, f);
                asian(i) = AsianOptionPrices(s0, v(i), T, r, sig, M2, f);
            end
        case 3
            for i=1:size(v, 2)
                [call(i), put(i)] = optionpricing(s0, K, T, v(i), sig, M, f);
                asian(i) = AsianOptionPrices(s0, K, T, v(i), sig, M2, f);
            end
        case 4
            for i=1:size(v, 2)
                [call(i), put(i)] = optionpricing(s0, K, T, r, v(i), M, f);
                asian(i) = AsianOptionPrices(s0, K, T, r, v(i), M2, f);
            end
    end
    plotting(call, put, asian, f, s, v);
end

function [call, put] = optionpricing(s0, K, T, r, sig, M, f)
    if(f==1)
        [u, d] = set1(T, M, sig);
    else
        [u, d] = set2(T, M, sig, r);
    end
    p  = (u-exp(r*(T/M)))/(u-d);
    %no arbitrage condition %
    if p<0||p>1
        fprintf('arbitrage detected function exiting ....%d, %d, %d, %d\n', s0, K, r, sig);
        call=0; 
        put=0;
        return;
    end
    s = finalStockPrices(s0, u, d, M);
    v = finalValuesCall(s, K, M);
    A = zeros(1, M+1);
    for i=0:M
        A(i+1) = v(i+1)*(p^i)*nchoosek(M, i);
        A(i+1) = A(i+1)*((1-p)^(M-i));
    end
    call = sum(A)*exp(-r*T);
    v = finalValuesPut(s, K, M);
    for i=0:M
        A(i+1) = v(i+1)*(p^i)*nchoosek(M, i);
        A(i+1) = A(i+1)*((1-p)^(M-i));
    end
    put  = sum(A)*exp(-r*T);  
end

function [asian] = AsianOptionPrices(s0, K, T, r, sig, M, f)
    if(f==1)
        [u, d] = set1(T, M, sig);
    else
        [u, d] = set2(T, M, sig, r);
    end
    p  = (u-exp(r*(T/M)))/(u-d);
    %no arbitrage condition %
    if p<0||p>1
        fprintf('arbitrage detected function exiting ....%d, %d, %d, %d\n', s0, K, r, sig);
        asian=0; 
        return;
    end
    v = finalValuesAsian(s0, K, u, d, M);
    for j=M:-1:1
        for i = 1:2:2^j-1
            k = floor((i+1)/2);
            v(1, k) = (p*v(1, i+1)+(1-p)*v(1, i))*(exp(-r*(T/M)));
        end
    end
    asian = v(1);
end
function [u, d] = set1(T, M, sig)
   u = exp(sig*sqrt(T/M));
   d = exp(-sig*sqrt(T/M));
end

function [u, d] = set2(T, M, sig, r)
   u = exp(sig*sqrt(T/M)+(r-0.5*sig*sig)*(T/M));
   d = exp(-sig*sqrt(T/M)+(r-0.5*sig*sig)*(T/M));
end

function v = finalValuesCall(s, K, M)
     v = zeros(1, M+1);
     for i=0:M
        v(1, i+1) = max(s(1, i+1)-K, 0);
     end
end
function v = finalValuesPut(s, K, M)
     v = zeros(1, M+1);
     for i=0:M
        v(1, i+1) = max(K-s(1, i+1), 0);
     end
end

function s = finalStockPrices(s0, u, d, M)
    s = zeros(1, M+1);
    for j=0:M
        s(1, j+1) = s0*(u^j)*(d^(M-j));
    end
    
end
function v = finalValuesAsian(s0, K, u, d, M)
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
            v2(1, 2*j-1) = s2(1, 2*j-1)+v(1, j);
            v2(1, 2*j) = s2(1, 2*j)+v(1, j);
        end
        for j=1:2^i
            v(1, j) = v2(1, j);
            s(1, j) = s2(1, j);
        end
    end
    for i=1:2^M
        v(1, i) = max([v(1,i)/(M+1)-K 0]);
    end
end
function varyM(s0, T, r, sig, M_v, M_v2, K_v2)
    for f=1:2
        for j=1:size(K_v2, 2)
            call = zeros(1, size(M_v, 2));
            put = zeros(1, size(M_v, 2));
            for k=1:size(M_v, 2)
                [call(k), put(k)] = optionpricing(s0, K_v2(j), T, r, sig, M_v(k), f);
            end
            asian = zeros(1, size(M_v2, 2));
            for k=1:size(M_v2, 2)
                asian(k) = AsianOptionPrices(s0, K_v2(j), T, r, sig, M_v2(k), f);
            end
            figure();
            plot(M_v, call);
            title(['Call Prices on varying parameter: M for set: ', num2str(f), ' and K = ', num2str(K_v2(j))]);
            xlabel('M');
            ylabel('Call Option Price');
            figure();
            plot(M_v, put);
            title(['Put Prices on varying parameter: M for set: ', num2str(f) , ' and K = ', num2str(K_v2(j))]);
            xlabel('M');
            ylabel('Put Option Price');
            figure();
            plot(M_v2, asian);
            title(['Asian Option Prices on varying parameter: M for set: ', num2str(f), ' and K = ', num2str(K_v2(j))]);
            xlabel('M');
            ylabel('Asian Option Price');
        end
    end
end
function plotting (call, put, asian, f, var, v)
    figure();
    hold on;
    plot(v, call);
    plot(v, put);
    title(['Option Prices on varying parameter: ', var, ' for set: ', num2str(f)]);
    xlabel(var);
    ylabel('Option Price');
    legend('Call', 'Put');
    hold off;
    figure();
    plot(v, asian);
    title(['Asian Option Prices on varying parameter: ', var, ' for set: ', num2str(f)]);
    xlabel(var);
    ylabel('Asian Option Price');
end

