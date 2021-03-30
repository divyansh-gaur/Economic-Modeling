
% This model simulates the growth of an economy based on the Swan-Solow model of economic growth with and without productivity shocks
 
%% Parameters
t = 250; % Number of years modelled in the simulation
n = 0.02; % The Population Grows at the rate of n per year 
s = 0.29; % Proportion of Income Saved or Invested
cons = 1-s; % Proportion of Income Consumed
b = 4.7; % Time discount factor
sigma = 3.0; % The Inverse of the elasticity of substitution of consumption
a = 0.4; % This is the elasticity of output with respect to capital (0<a<1)
d = 0.073;  % Stock of Capital Depreciates at the rate of d every year
g = 0.01; % Technology Grows at the rate of g per year
rho = 0.88;  % Autoregressive parameter in Productivity shock
sigma_e = 0.01; % Standard Deviation of the innovation in the autoregressive process for the technology shock

% Steady State Parameters
K_ss = (s*b/((1+n)*(1+g)-(1-d)))^(1/(1-a)); % Steady State stock of capital
Y_ss = b*(K_ss^a); % Output in Steady State
C_ss = cons*Y_ss; % Consumption in Steady State
I_ss = s*Y_ss; % Investment (Savings) in Steady State
marg_productivity_ss = a*b*K_ss^(a-1); % Marginal Productivity

%% Initial Values
Y0 = 21427675; % Total Output in Year 0
K0 = 56215.312; % Initial Stock of Capital
labor_force = 253.069; % Initial Labor Force of the country
A0 = Y0/( K0^(a)*labor_force^(1-a));

%% Simulation

% Initializing the matrices
A = zeros();
N = zeros();
K = zeros();

K(1) = K0/labor_force;
N(1) = labor_force;
A(1) = A0; 

% Recursive Algorithm that models the economy
for j = 1:t-1
    N(j+1,1) = (1+n)*N(j,1);
    A(j+1,1) = (1+g)*A(j,1);
    K(j+1,1) = (1-d)*K(j,1)/(1+n) + s*A(j,1)*K(j,1)^(a)/(1+n);
end

%% Variable Computation

Y = A.*(K.^(a)).*(N.^(1-a)); % Total demand/output in the economy
G = A.*K.^(a)*1000; % Effective output per worker in the economy

K = K.*N; % Capital Stock in the economy
cost_cap = a*A.*(K.^(a-1)).*(N.^(1-a)); % The cost of capital/cost of borrowing

real_wages = (1-a)*A.*(K.^(a)).*(N.^(-a)); % Real wages earned by the workers

%% Plotting 

figure(1)
subplot(3, 2, 1)
plot((0:(t-1)),Y, 'r')
hold on
plot((0:(t-1)), cons.*Y, 'g')
plot((0:(t-1)), s.*Y, 'b')
title('Real GDP and Household Breakdown');
xlabel('Time from Year 0');
ylabel('Dollars (Billions)');
legend('Real GDP','Consumption (Adjusted for Inflation)','Investment (Adjusted for Inflation)', 'Location', 'Best')
hold off

subplot(3, 2, 2)
semilogy((0:(t-1)),Y, 'b')
title('Semilog Real GDP');
xlabel('Time from Year 0');
ylabel('Dollars (Billions)');

subplot(3, 2, 3)
plot((0:(t-1)),G, 'm')
title('Output per worker (Adjusted for Inflation)');
xlabel('Time from Year 0');
ylabel('Inflation Adjusted Output');

subplot(3, 2, 4)
plot((0:(t-1)),cost_cap, 'g')
title('Cost of Capital');
xlabel('Time from Year 0');
ylabel('Real Cost of Capital');

subplot(3, 2, 5)
plot((0:(t-1)),K, 'k')
title('Capital-to-Labor Ratio in the long run');
xlabel('Time from Year 0');
ylabel('Capital-to-Labor Ratio');

subplot(3, 2, 6)
plot((0:(t-1)),real_wages, 'c')
title('Wages (Adjusted for Inflation) in the Long Run');
xlabel('Time from Year 0');
ylabel('Real Wages');

pause

%% Shock Simuation

% Initial Values
theta_0 = 1; % Gives the output to capital ratio in the economy
K0_shock = K_ss; 
Y0_shock = b*theta_0*(K0_shock^a);
C0_shock = (1-s)*Y0_shock; 
I0_shock = s*Y0_shock; 
marg_productivity0 = a*b*theta_0*K0_shock^(a-1);

% Initializing the matrices
theta = []; % Theta Vector
K_shock = []; % Capital Vector
Y_shock = []; % Output Vector
C_shock = []; % Consumption Vector
I_shock = []; % Invertment Vector
marg_productivity = []; % Marginal Productivity Vector
Y_grow = []; % Growth in Output Vector

K_shock = [K_shock; K0_shock];
Y_shock = [Y_shock; Y0_shock];
C_shock = [C_shock; C0_shock];
I_shock = [I_shock; I0_shock];
marg_productivity = [marg_productivity; marg_productivity0];
theta = [theta; theta_0];

% Introduction of productivity shock through change in technology
tech = randn(t+100,1)*sigma_e; 
tech = tech(101:t+100);

% Recursive Algorithm that models the economy
for p = 1:t
    
    theta_p = exp(rho*log(theta_0)+tech(p)); 
    theta_0 = theta_p; 
    theta = [theta; theta_p];
    
    K_p = (s/((1+n)*(1+g)))*b*theta_p*(K0_shock^a)+((1-d)/((1+n)*(1+g)))*K0_shock;
    K0_shock = K_p;
    K_shock = [K_shock; K_p];
    
    Y_p = b*theta_p*K_p^a;
    Y_shock = [Y_shock; Y_p];
    
    C_p = (1-s)*Y_p;
    C_shock = [C_shock; C_p];
    
    I_p = s*Y_p;
    I_shock = [I_shock; I_p];
    
    marg_productivity_p = a*b*theta_p*K_p^(a-1);
    marg_productivity = [marg_productivity; marg_productivity_p];
    
end

% Computation of Logarithmic Output
log_output = 1:t+1; 
log_output = log_output';
ln_Y = log(Y_shock)+log(1+g).*log_output;

% Computation of Output Growth
for q = 1:t
    Y_grow_t = ((Y_shock(q+1)-Y_shock(q))/Y_shock(q))*100;
    Y_grow = [Y_grow; Y_grow_t];
end
%% Plotting for Shocks

figure(2)
suptitle('Productivity Shocks')
subplot(3,2,1)
plot(Y_shock(1:t), 'k');
title('Real GDP (Output)');
xlabel('Time from Year 0');

subplot(3,2,2)
plot(ln_Y(1:100));
title('Output on a Logarithmic Scale');
xlabel('Time from Year 0');

subplot(3,2,3)
plot(Y_grow(1:t), 'b');
title('Growth rate of Output');
xlabel('Time from Year 0');

subplot(3,2,4)
plot(C_shock(1:t), 'g');
hold on
plot(I_shock(1:t), 'y');
title('Total Consumption & Investment (Adjusted for Inflation');
xlabel('Time from Year 0');
legend('Consumption', 'Investment');
hold off

subplot(3,2,5)
plot(K_shock(1:t), 'r');
title('Stock of capital');
xlabel('Time from Year 0');

subplot(3,2,6)
plot(marg_productivity(1:t), 'c');
title('Real interest rate');
xlabel('Time from Year 0');
