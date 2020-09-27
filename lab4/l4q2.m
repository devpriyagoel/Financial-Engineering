clc;
close all;
warning('off','all')

m = [0.1 0.2 0.15];
C = [ 0.005 -0.010 0.004; -0.010 0.040 -0.002; 0.004 -0.002 0.023];
u = ones(1, size(m, 2));

eff_port_return = 0:0.0001:1;
asset_risk = sqrt(diag(C));

eff_port_risk = evaluate_min_portfolio_risk_with_given_return(eff_port_return, C, u, m);
mvp = minimum_variance_portfolio(m, u, C);
eff_mvp_port_return = mvp:.001:1;
eff_mvp_port_risk = evaluate_min_portfolio_risk_with_given_return(eff_mvp_port_return, C, u, m);

% For 1 and 2 
C1 = [ 0.005 -0.010; -0.010 0.040];
m1 = [0.1 0.2];
u1 = ones(1, size(m1, 2));
eff_port_risk1 = evaluate_min_portfolio_risk_with_given_return(eff_port_return, C1, u1, m1);
mvp1 = minimum_variance_portfolio(m1, u1, C1);
eff_mvp_port_return1 = mvp1:.001:1;
eff_mvp_port_risk1 = evaluate_min_portfolio_risk_with_given_return(eff_mvp_port_return1, C1, u1, m1);

% For 2 and 3
C2 = [0.040 -0.002;-0.002 0.023];
m2 = [0.2 0.15];
u2 = ones(1, size(m2, 2));
eff_port_risk2 = evaluate_min_portfolio_risk_with_given_return(eff_port_return, C2, u2, m2);
mvp2 = minimum_variance_portfolio(m2, u2, C2);
eff_mvp_port_return2 = mvp2:.001:1;
eff_mvp_port_risk2 = evaluate_min_portfolio_risk_with_given_return(eff_mvp_port_return2, C2, u2, m2);

% For 1 and 3 
C3 = [ 0.005 0.004; 0.004 0.023];
m3 = [0.1 0.15];
u3 = ones(1, size(m3, 2));
eff_port_risk3 = evaluate_min_portfolio_risk_with_given_return(eff_port_return, C3, u3, m3);
mvp3 = minimum_variance_portfolio(m3, u3, C3);
eff_mvp_port_return3 = mvp3:.001:1;
eff_mvp_port_risk3 = evaluate_min_portfolio_risk_with_given_return(eff_mvp_port_return3, C3, u3, m3);

l = 100000;
r1 = unifrnd(0.0, 10.0, l, 1);
r2 = unifrnd(0.0, 10.0, l, 1);
r3 = unifrnd(0.0, 10.0, l, 1);

reg_risk = zeros(1, l);
reg_return = zeros(1, l);
for i=1:l
    w = [r1(i) r2(i) r3(i)]/(r1(i)+r2(i)+r3(i));
    reg_risk(i) = portfolio_risk(w, C);
    reg_return(i) = portfolio_return(m, w);
end
figure();
hold on;
grid on;
plot(eff_port_risk, eff_port_return, 'lineWidth',2.5,'color','c');
plot(eff_mvp_port_risk, eff_mvp_port_return,'lineWidth',1.5,'color','b');
plot(eff_port_risk1, eff_port_return, 'lineWidth',2.5,'color','y');
plot(eff_mvp_port_risk1, eff_mvp_port_return1,'lineWidth',1.5,'color','g');
plot(eff_port_risk2, eff_port_return, 'lineWidth',2.5,'color','#D95319');
plot(eff_mvp_port_risk2, eff_mvp_port_return2,'lineWidth',1.5,'color','r');
plot(eff_port_risk3, eff_port_return, 'lineWidth',2.5,'color','#7E2F8E');
plot(eff_mvp_port_risk3, eff_mvp_port_return3,'lineWidth',1.5,'color','k');
plot(asset_risk, m, '*','color','#A2142F');
plot(reg_risk, reg_return, '.', 'color', '#7f5b60');

xlim([0, 0.5]);
title('Minimum Variance Curve, Efficient Frontier and Feasible Region');
ylabel('Expected Return \mu');
xlabel('Standard Deviation \delta');
legend('Minimum Variance Curve','Efficient Frontier', 'Minimum Variance Curve for 1&2','Efficient Frontier for 1&2','Minimum Variance Curve for 2&3','Efficient Frontier for 2&3','Minimum Variance Curve for 1&3','Efficient Frontier for 1&3','Individual Assets', 'Feasible Region');

k=1;
for i= 1:size(eff_port_return, 2)
    w = minimum_risk_portfolio_with_given_return(eff_port_return(i), C, u, m);
    if(w(1)>=0&&w(2)>=0&&w(3)>=0)
        X(k) = w(1);
        Y(k) = w(2);
        k = k+1;
    end
    w = minimum_risk_portfolio_with_given_return(eff_port_return(i), C1, u1, m1);
    if(w(1)>=0&&w(2)>=0)
        X(k) = w(1);
        Y(k) = w(2);
        k = k+1;
    end
    w = minimum_risk_portfolio_with_given_return(eff_port_return(i), C2, u2, m2);
    if(w(1)>=0&&w(2)>=0)
        X(k) = 1-w(1)-w(2);
        Y(k) = w(1);
        k = k+1;
    end
    w = minimum_risk_portfolio_with_given_return(eff_port_return(i), C3, u3, m3);
    if(w(1)>=0&&w(2)>=0)
        X(k) = w(1);
        Y(k) = 1-w(1)-w(2);
        k = k+1;
    end
end
figure();
plot(X, Y, '.');
title('Weights w2 vs w1 corresponding to minimum variance curve eq: w1+w2+w3=1');
ylabel('w2');
xlabel('w1');

fprintf('Equation satisfied by weights: w1+w2+w3=1\n');
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