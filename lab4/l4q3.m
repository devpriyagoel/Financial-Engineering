clc;
close all;
warning('off','all')

datafile='./Data.csv';
data=csvread(datafile,1,1);
m=(mean(data));
C=cov(data);
u = ones(1, size(m, 2));


%part a%
disp('Part A');
eff_port_return = 0:0.1:250;
eff_port_risk = evaluate_min_portfolio_risk_with_given_return(eff_port_return, C, u, m);
asset_risk = sqrt(diag(C));
mvp = minimum_variance_portfolio(m, u, C);
eff_mvp_port_return = mvp:.1:250;
eff_mvp_port_risk = evaluate_min_portfolio_risk_with_given_return(eff_mvp_port_return, C, u, m);
figure();
hold on;
grid on;
plot(eff_port_risk, eff_port_return, 'lineWidth',2.5,'color','b');
plot(eff_mvp_port_risk, eff_mvp_port_return,'lineWidth',1.5,'color','m');
%plot(asset_risk, m, '*','color','k');
xlim([0.5, 1.8]);
ylim([0, 250]);
title('Markowitz efficient frontier for different assets');
ylabel('Expected Return \mu');
xlabel('Standard Deviation \delta');
legend('Minimum Variance Curve','Efficient Frontier');
hold off;
%part b%
fprintf('\n');
disp('Part B');
murf = 0.05;
wcap=(((m-murf).*u)/C)/((((m-murf).*u)/C)*u');
rcap = portfolio_return(m, wcap);
sigcap = portfolio_risk(wcap, C);

fprintf('For 5 %% risk-free return:\n');    
fprintf('Weights of Market portfolio: %6.4f  %6.4f  %6.4f\n',wcap(1,1),wcap(1,2),wcap(1,3));
fprintf('Return on Market portfolio: %6.4f\n',rcap);
fprintf('Risk on Market portfolio: %6.4f\n',sigcap);

%part c%
fprintf('\n');
disp('Part C');
x=0.5:0.1:1.8;
y=murf+((rcap-murf)/sigcap)*x;
figure();
hold on;
grid on;
plot(x,y,'color','r', 'linewidth', 1.5);
plot(eff_mvp_port_risk, eff_mvp_port_return,'lineWidth',1,'color','g');
title('Markowitz Efficient Frontier & Capital Market Line');
xlim([0.5, 1.8]);
ylim([0, 250]);
xlabel('Risk (\sigma)');
ylabel('Return (\mu)');
plot(sigcap,rcap,'*');
legend('CML','Markowitz Efficient Frontier','Market Portfolio');
syms x;
y=murf+vpa(((rcap-murf)/sigcap))*x;
disp('Capital Market Line:');
disp(y);
%part d%
fprintf('\n');
disp('Part D');
beta=-2:0.1:2;
muv=murf+(rcap-murf).*beta;
figure();
plot(beta,muv);
grid on
title('Security Market Line')
xlabel('Beta (\beta)');
ylabel('Return (\mu)');
syms x;
y = murf+vpa((rcap-murf)).*x;
disp('Security Market Line:');
disp(y);
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

