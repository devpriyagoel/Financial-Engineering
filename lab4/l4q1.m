clc;
close all;
warning('off','all')

m = [0.1 0.2 0.15];
C = [ 0.005 -0.010 0.004; -0.010 0.040 -0.002; 0.004 -0.002 0.023];
u = ones(1, size(m, 2));


%part a%
disp('Part A');
eff_port_return = 0:0.001:1;
eff_port_risk = evaluate_min_portfolio_risk_with_given_return(eff_port_return, C, u, m);
asset_risk = sqrt(diag(C));
mvp = minimum_variance_portfolio(m, u, C);
eff_mvp_port_return = mvp:.001:1;
eff_mvp_port_risk = evaluate_min_portfolio_risk_with_given_return(eff_mvp_port_return, C, u, m);
figure();
hold on;
grid on;
plot(eff_port_risk, eff_port_return, 'lineWidth',2.5,'color','b');
plot(eff_mvp_port_risk, eff_mvp_port_return,'lineWidth',1.5,'color','m');
plot(asset_risk, m, '*','color','k');

title('Markowitz efficient frontier for different assets');
ylabel('Expected Return \mu');
xlabel('Standard Deviation \delta');
legend('Minimum Variance Curve','Efficient Frontier', 'Individual Assets');

%part b%
fprintf('\n');
disp('Part B');
fprintf('   Return   Risk\t  W1\t\tW2\t\tW3\n');
for i=0:0.1:1 
    w = minimum_risk_portfolio_with_given_return(i, C, u, m);
    eff_port_risk = portfolio_risk(w, C);
    fprintf('%8.2f %8.4f %8.4f %8.4f %8.4f\n',i,eff_port_risk,w(1,1),w(1,2),w(1,3));
end

%part c%
fprintf('\n');
disp('Part C');
eff_port_risk = 0.15;

[h1, h2] = evaluate_return(C, u, m,eff_port_risk);

w = minimum_risk_portfolio_with_given_return(h1, C, u, m);
fprintf('Maximum Return Portfolio for 15 %% risk: \n');
fprintf('Return: %6.4f and Weights: %6.4f  %6.4f  %6.4f\n',h1,w(1,1),w(1,2),w(1,3));

w = minimum_risk_portfolio_with_given_return(h2, C, u, m);
fprintf('Minimum Return Portfolio for 15 %% risk:\n');
fprintf('Return: %6.4f and Weights: %6.4f  %6.4f  %6.4f\n',h2,w(1,1),w(1,2),w(1,3));

%part d%
fprintf('\n');
disp('Part D');
h = 0.18;
w = minimum_risk_portfolio_with_given_return(h, C, u, m);
sig = portfolio_risk(w, C);
fprintf('Minimum Risk Portfolio for 18 %% return:\n');    
fprintf('Risk: %6.4f and Weights: %6.4f  %6.4f  %6.4f\n',sig,w(1,1),w(1,2),w(1,3));

%part e%
fprintf('\n');
disp('Part E');
murf = 0.1;
wcap=(((m-murf).*u)/C)/((((m-murf).*u)/C)*u');
rcap = portfolio_return(m, wcap);
sigcap = portfolio_risk(wcap, C);

fprintf('For 10 %% risk-free return:\n');    
fprintf('Weights of Market portfolio: %6.4f  %6.4f  %6.4f\n',wcap(1,1),wcap(1,2),wcap(1,3));
fprintf('Return on Market portfolio: %6.4f\n',rcap);
fprintf('Risk on Market portfolio: %6.4f\n',sigcap);

x=0:0.1:1.8;
y=murf+((rcap-murf)/sigcap)*x;
figure();
hold on;
grid on;
plot(eff_mvp_port_risk, eff_mvp_port_return, 'lineWidth',1,'color','g');
plot(x,y,'color','r', 'linewidth', 1.5);
title('Capital Market Line');
xlabel('Risk (\sigma)');
ylabel('Return (\mu)');
plot(sigcap,rcap,'*', 'color', 'k');
legend('Markowitz Efficient Frontier','CML','Market Portfolio');
syms x;
y=murf+vpa(((rcap-murf)/sigcap))*x;
disp('Capital Market Line:');
disp(y);

%Part f
fprintf('\n');
disp('Part F');
x = 0.1;
y=murf+((rcap-murf)/sigcap)*x;
riskfree=(y-rcap)/(murf-rcap);
risky=(1-riskfree).*wcap;

fprintf('Portfolio for 10 %% risk:\n');  
fprintf('Risk-free Asset Weight: %6.4f\n',riskfree);
fprintf('Risky Assets Weights: %6.4f  %6.4f  %6.4f\n',risky(1,1),risky(1,2),risky(1,3));

x=0.25;
y=murf+((rcap-murf)/sigcap)*x;
riskfree=(y-rcap)/(murf-rcap);
risky=(1-riskfree).*wcap;

fprintf('Portfolio for 25 %% risk:\n');  
fprintf('Risk-free Asset Weight: %6.4f\n',riskfree);
fprintf('Risky Assets Weights: %6.4f  %6.4f  %6.4f\n',risky(1,1),risky(1,2),risky(1,3));


function mu = portfolio_return(m, w)
    mu = m*w';
end
function sig = portfolio_risk(w, C)
    sig = sqrt(w*C*w');
end
function ret = minimum_variance_portfolio(m, u, C)
    w = u/C;
    w = w/(w*u');
    ret = portfolio_return(m, w);
end
function w = minimum_risk_portfolio_with_given_return(port_return, C, u, m)
    uCinv = u/C;
    mCinv = m/C;
    a = uCinv*u';
    b = uCinv*m';
    c = mCinv*u';
    d = mCinv*m';
    det1 = evaluate_det(a, b, c, d);
    det2 = evaluate_det(1, b, port_return, d);
    det3 = evaluate_det(a, 1, c, port_return);
    w = (det2*uCinv)/det1 + (det3*mCinv)/det1;
end
function port_risk = evaluate_min_portfolio_risk_with_given_return(port_return, C, u, m)
    port_risk = zeros(1, size(port_return, 2));
    for i=1:size(port_return, 2)
        w = minimum_risk_portfolio_with_given_return(port_return(i), C, u, m);
        port_risk(i) = portfolio_risk(w, C);
    end
end

function det = evaluate_det(a, b, c, d)
    det = a*d-b*c;
end
function [r1, r2] = evaluate_return(C, u, m, sg)
    uCinv = u/C;
    mCinv = m/C;
    a = uCinv*u';
    b = uCinv*m';
    c = mCinv*u';
    d = mCinv*m';
    det = evaluate_det(a, b, c, d);
    r1=(2*b+sqrt(4*b*b-4*a*(d-det*sg*sg)))/(2*a);
    r2=(2*b-sqrt(4*b*b-4*a*(d-det*sg*sg)))/(2*a);
end

