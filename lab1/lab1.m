clc;
clear all;
close all;
warning('off','all')

s0 = 100;
K = 105;
T = 5;
r = 0.05;
q = 0.3;

%part 1%
M = [1 5 10 20 50 100 200 400];
[call, put] = initialOptionPrices(M, s0, K, T, r, q);
fprintf('QUESTION 1: Initial Option Prices for given number of subintervals\n');
fprintf('  M\t Call prices \t Put Prices\n');
for i=1:size(M, 2)
    fprintf('%3d\t %8.5f\t %8.5f\n', M(i), call(i), put(i));
end

%part 2%
fprintf('QUESTION 2: Initial Option Prices Vs number of subintervals plots\n');
M = 1:100;
[call, put] = initialOptionPrices(M, s0, K, T, r, q);
figure();
plot(M, call);
title('Call Option for step size 1');
xlabel('M');
ylabel('Option Price');
figure();
plot(M, put);
title('Put Option for step size 1');
xlabel('M');
ylabel('Option Price');
M=1:5:200;
[call, put] = initialOptionPrices(M, s0, K, T, r, q);
figure();
plot(M, call);
title('Call Option for step size 5');
xlabel('M');
ylabel('Option Price');
figure();
plot(M, put);
title('Put Option for step size 5');
xlabel('M');
ylabel('Option Price');

[call, put] = BSMcall(0, s0, q, T, K, r);
disp('Call Option Price Using BSM formula for continous Model');
disp(call);
disp('Put Option Price Using BSM formula for continous Model');
disp(put);
%part 3%
w = [0 0.5 1 1.5 3 4.5];
M = 20;
l = (M*w)/T;
[call, put, stock] = intermediatePrices(s0, K, T, r, q, M, w, l);
fprintf('QUESTION 3: Option prices at given time values for the case M = 20: \n');
fprintf('\nOption prices: \n');
for i=1:size(w, 2)
    fprintf('\n t = %5.3f\n', w(i));
    fprintf('Stock Price\t Call Price\t Put Price\n');
    for j=0:l(i)
        fprintf('%8.5f\t %8.5f\t %8.5f\n', stock(i, j+1), call(i, j+1), put(i, j+1));
    end
end
fprintf('\n');


function [call, put] = optionpricing(s0, K, T, r, q, M)
    u = exp(q*sqrt(T/M)+(r-0.5*q*q)*(T/M));
    d = exp(-q*sqrt(T/M)+(r-0.5*q*q)*(T/M));
    p  = (u-exp(r*(T/M)))/(u-d);
    %no arbitrage condition %
    if p<0||p>1
        fprintf('arbitrage detected function exiting ....\n');
        return;
    end
    s = finalStockPrices(s0, u, d, M);
    v = finalValues(s, K, M);
    A = zeros(1, M+1);
    for i=0:M
        A(i+1) = v(i+1)*(p^i)*nchoosek(M, i);
        A(i+1) = A(i+1)*((1-p)^(M-i));
    end
    call = sum(A)*exp(-r*T);
    put  = call + K*exp(-r*T)-s0;  %put call parity%
end

function [call, put, stock] = intermediatePrices(s0, K, T, r, q, M, w, l)
    u = exp(q*sqrt(T/M)+(r-0.5*q*q)*(T/M));
    d = exp(-q*sqrt(T/M)+(r-0.5*q*q)*(T/M));
    p  = (u-exp(r*(T/M)))/(u-d);
    %no arbitrage condition %
    if p<0||p>1
        fprintf('arbitrage detected function exiting ....\n');
        return;
    end
    [s, v] = portfolioPriceTree(s0, K, T, r, u, d, p, M);
    call = zeros(size(w, 2), M+1);
    put = zeros(size(w, 2), M+1);
    stock = zeros(size(w, 2), M+1);
    for i = 1:size(w, 2)
        for j = 0:l(i)
            call(i, j+1) = v(l(i)+1, j+1);
            stock(i, j+1) = s(l(i)+1, j+1); 
            put(i, j+1) = call(i, j+1) + K*exp(-r*(T-w(i)))-stock(i, j+1); %put call parity%
        end
    end
end

function v = finalValues(s, K, M)
     v = zeros(1, M+1);
     for i=0:M
        v(1, i+1) = max(s(1, i+1)-K, 0);
     end
end
function s = finalStockPrices(s0, u, d, M)
    s = zeros(1, M+1);
    for j=0:M
        s(1, j+1) = s0*(u^j)*(d^(M-j));
    end
    
end
function [s, v] = portfolioPriceTree(s0, K, T, r, u, d, p, M)
    s = zeros(M+1, M+1);
    v = zeros(M+1, M+1);
    for i=0:M
        for j=0:i
            s(i+1, j+1) = s0*(u^j)*(d^(i-j));
        end
    end
    v(M+1, :) = finalValues(s(M+1, :), K, M);
    for i=M-1:-1:0
        for j=0:i
            v(i+1, j+1) = (p*v(i+2, j+1)+(1-p)*v(i+2, j+2))*(exp(-r*(T/M)));
        end
    end
end
function [call, put] = initialOptionPrices(M, s0, K, T, r, q)
    call = zeros(1, size(M, 2));
    put = zeros(1, size(M, 2));
    for i=1:size(M, 2)
        [call(i), put(i)] = optionpricing(s0, K, T, r, q, M(i));
    end
end


function [call, put] = BSMcall(t, x, q, T, K, r)
    call = x*normcdf((log(x/K)+r*(T-t)+0.5*q*q*(T-t))/(q*sqrt(T-t)))-exp(-r*(T-t))*K*normcdf((log(x/K)+r*(T-t)-0.5*q*q*(T-t))/(q*sqrt(T-t)));
    put = call + K*exp(-r*T)-x;
end