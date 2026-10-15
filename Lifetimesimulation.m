% lifetime_simulation_4panel.m
% -------------------------------------------------------------------------
% LIFETIME SIMULATION (FIGURE 8 & TABLE 2)
% Biological Purpose: Simulates human ageing from 30 to 100 years old by 
% applying a linear age drift to baseline amyloid production. It tracks all 
% four state variables and precisely identifies the biological tipping point.
% -------------------------------------------------------------------------

clear; clc; close all;

%% ========================================================================
% 1. PARAMETER INITIALIZATION
% ========================================================================
% Neuronal death hazard parameters
sr = 0.20; d0 = 2e-4; da = 0.002; dtau = 0.05;
Vtox = 0.12; Ktox = 0.4;

% Microglial & Switch parameters
rm = 1.0; beta = 4.0; n = 4; sigmaStar = 0.75;
Ksw = sigmaStar/(1-sigmaStar);          % Half-saturation constant (K = 3)
gamma = 0.5;

% Amyloid-beta parameters
pstim = 2.5; Amax = 20; cclear = 1.0; Vc = 10; theta = 0.3;

% Tau dynamics parameters
ktau = 0.1; tau0 = 0.05; btau = 0.5; gtau = 2.0; dt = 0.1;

% Age Drift Parameters (Human Ageing Engine)
% Amyloid production increases linearly: p(t) = 2.0 + 0.1*t
p0 = 2.0;       % Baseline A-beta production at age 30
alpha = 0.1;    % Linear annual drift coefficient

% Initial Conditions at Age 30 (Healthy Baseline)
% [s=1 (100% neurons), m=0 (no toxic microglia), a=0.2 (low amyloid), tau=0.05]
x0 = [1; 0; 0.2; 0.05];

% High-precision solver options
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.05);

%% ========================================================================
% 2. ODE SYSTEM SETUP & INTEGRATION
% ========================================================================
% Define the full 4-ODE mathematical model
rhs = @(t,x) lifetime(t,x,sr,d0,da,dtau,Vtox,Ktox,rm,beta,n,Ksw,gamma, ...
    pstim,Amax,cclear,Vc,theta,ktau,tau0,btau,gtau,dt,p0,alpha);

% Define Control 3 (Vtox = 0) to show survival if microglia couldn't kill
rhs0 = @(t,x) lifetime(t,x,sr,d0,da,dtau,0,Ktox,rm,beta,n,Ksw,gamma, ...
    pstim,Amax,cclear,Vc,theta,ktau,tau0,btau,gtau,dt,p0,alpha);

% Integrate both systems over 70 years (from t=0 to t=70)
[t,X] = ode45(rhs, [0 70], x0, opts);
[t0,X0] = ode45(rhs0, [0 70], x0, opts);

% Translate simulation time (t) into human age (t=0 means age 30)
age = t + 30;
age0 = t0 + 30;

% Extract the state variables for the full model
s = X(:,1); m = X(:,2); a = X(:,3); tau = X(:,4);

%% ========================================================================
% 3. IDENTIFYING THE BIOLOGICAL SWITCH
% ========================================================================
% Mathematically define the switch as the point where the toxic microglial 
% fraction is increasing at its absolute fastest rate (maximum dm/dt).
H = a.^n ./ (Ksw^n + a.^n);
mdot = rm .* ((1-m).*(beta*H + gamma*m) - m);
[~,ix] = max(mdot); % Find the array index of the maximum rate of change
ageSw = age(ix);    % Translate that index into the precise human age

% Print the switch age and the production rate at that exact moment
fprintf('microglial switch (max m-dot) at age %.1f, p = %.2f\n', ageSw, p0 + alpha*(ageSw-30));

%% ========================================================================
% 4. GENERATING TABLE 2 (CONSOLE OUTPUT)
% ========================================================================
% This loop interpolates the ODE data at specific decade checkpoints to 
% print the values required for "Table 2: Lifetime trajectory" in the report.
chk = [30 40 50 60 70 75 80 85 90 100];
fprintf('\n%-6s %8s %8s %8s %8s\n', 'age', 's', 'm', 'a', 'tau');
for k = 1:numel(chk)
    tk = chk(k) - 30; % Convert target age back to simulation time 't'
    fprintf('%-6.0f %8.3f %8.3f %8.3f %8.3f\n', chk(k), ...
        interp1(t,s,tk), interp1(t,m,tk), interp1(t,a,tk), interp1(t,tau,tk));
end

%% ========================================================================
% 5. GENERATING FIGURE 8 (4-PANEL LIFETIME TRAJECTORY)
% ========================================================================
figure('Color','w','Position',[80 60 860 640]);

% Setup labels, colors, and data for the 2x2 subplot layout
yl = {'surviving neurons s', 'toxic microglia m', 'A\beta load a', 'tangle burden \tau'};
cols = {[0 0 0], [0.85 0.45 0], [0 0.35 0.85], [0.15 0.55 0.25]};
series = {s, m, a, tau};

for ip = 1:4
    subplot(2,2,ip); hold on; box on;
    set(gca, 'FontSize', 11, 'Layer', 'top');
    
    % Special additions for Panel 1 (Surviving Neurons)
    if ip == 1
        % Plot the Vtox=0 control run (dotted black line)
        plot(age0, X0(:,1), 'k:', 'LineWidth', 1.4);
        % Plot the mathematical resistant floor (s_r = 0.20)
        yline(sr, '-.', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2);
    end
    
    % Plot the main state variable trajectory
    plot(age, series{ip}, 'Color', cols{ip}, 'LineWidth', 1.8);
    
    % Add the red vertical dashed line indicating the microglial switch
    xline(ageSw, 'r--', 'LineWidth', 1.2);
    
    % Special addition for Panel 3 (Amyloid Load)
    if ip == 3
        % Mark the K=3 threshold where microglial clearance begins to fail
        yline(Ksw, 'r:', 'LineWidth', 1.2);
    end
    
    % Formatting
    xlim([30 100]);
    xlabel('age (yr)');
    ylabel(yl{ip});
end

% Add legend specifically to the first subplot (Neurons)
subplot(2,2,1);
legend({'V_{tox}=0', 's_r = 0.20', 'full model', 'switch'}, 'Location', 'southwest');

% Main figure title
sgtitle('Lifetime simulation, p(t) = 2.0 + 0.1 t yr^{-1}');

%% ========================================================================
% 6. HELPER FUNCTION: THE FULL MATHEMATICAL ODE SYSTEM
% ========================================================================
% Evaluates the 4 differential equations at every time step 't'
function dx = lifetime(t, x, sr, d0, da, dtau, Vtox, Ktox, ...
        rm, beta, n, K, gamma, pstim, Amax, cclear, Vc, theta, ...
        ktau, tau0, btau, gtau, dt, p0, alpha)
    
    % Unpack state variables
    s = x(1); m = x(2); a = x(3); tau = x(4);
    
    % Upstream functions
    sig = a/(1+a);
    H = a^n/(K^n + a^n);
    
    % Apply the linear age drift to baseline amyloid production
    p = p0 + alpha*t;
    
    % Calculate total per-capita neuronal death hazard
    Phi = d0 + da*sig + dtau*tau^2/(1+tau^2) + Vtox*m^2/(Ktox^2 + m^2);
    
    % Return the column vector [ds/dt; dm/dt; da/dt; dtau/dt]
    dx = [ -(s-sr)*Phi;
           rm*((1-m)*(beta*H + gamma*m) - m);
           (p + pstim*m)*(1 - a/Amax) - cclear*a - Vc*((1-m)+theta*m)*sig;
           ktau*s*(tau0 + btau*sig + gtau*m) - dt*tau ];
end