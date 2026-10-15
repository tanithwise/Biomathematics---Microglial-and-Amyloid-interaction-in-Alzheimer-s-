% control3_neuronal_survival_age_drift.m
% -------------------------------------------------------------------------
% CONTROL 3: MICROGLIAL TOXICITY ONLY & LIFETIME SIMULATION
% Biological Purpose: Simulates the chronological progression of the disease 
% from age 30 to 100. It isolates the direct contribution of microglial 
% killing (Vtox) to neuronal collapse while keeping the fold unchanged.
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

% Age Drift Parameters (Human Ageing)
% Amyloid production increases linearly: p(t) = p0 + alpha*t
p0 = 2.0;       % Baseline production at age 30
alpha = 0.1;    % Linear annual drift coefficient

% Initial Conditions at Age 30 (s=1, m=0, a=0.2, tau=0.05)
x0 = [1; 0; 0.2; 0.05];

% High-precision solver options for ode45
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.05);

%% ========================================================================
% 2. SIMULATION CASES DEFINITION
% ========================================================================
% We test three specific scenarios to isolate the killing mechanism.
runs = { ...
    % 1. The full biological model with all toxic feedback loops active
    'full model',            theta, pstim, Vtox; ...
    % 2. Control 3: The switch happens, but microglia cannot directly kill neurons
    'V_{tox} = 0',           theta, pstim, 0; ...
    % 3. Control 2: The toxic microglial phenotype is completely disabled
    'toxic switch removed',  1,     0,     0};

%% ========================================================================
% 3. FIGURE SETUP
% ========================================================================
figure('Color','w','Position',[120 80 560 420]);
hold on; box on;
set(gca, 'FontSize', 11, 'Layer', 'top');

styles = {'k-', 'k:', 'Color'};
styleColor = [0.15 0.55 0.25];
ageSw = nan; % Variable to store the exact age of the microglial switch

%% ========================================================================
% 4. LIFETIME INTEGRATION (AGE 30 TO 100)
% ========================================================================
for i = 1:3
    % Extract the specific parameters for this run
    th = runs{i,2}; ps = runs{i,3}; Vt = runs{i,4};
    
    % Define the 4-ODE system utilizing the lifetime helper function
    rhs = @(t,x) lifetime(t,x,sr,d0,da,dtau,Vt,Ktox,rm,beta,n,Ksw,gamma, ...
        ps,Amax,cclear,Vc,th,ktau,tau0,btau,gtau,dt,p0,alpha);
    
    % Integrate the system over 70 years (t=0 represents age 30)
    [t,X] = ode45(rhs, [0 70], x0, opts);
    
    % Convert simulation time (t) to human age
    age = t + 30;
    
    % Plotting and Data Extraction based on the current run
    if i == 1
        % Plot the full model (solid black)
        plot(age, X(:,1), 'k-', 'LineWidth', 1.8);
        
        % Mathematically define the Microglial Switch:
        % The switch occurs when the toxic fraction is increasing the fastest (max dm/dt).
        H = X(:,3).^n ./ (Ksw^n + X(:,3).^n);
        mdot = rm .* ((1-X(:,2)).*(beta*H + gamma*X(:,2)) - X(:,2));
        [~,ix] = max(mdot); % Find the index of the maximum rate of change
        ageSw = age(ix);    % Translate that index into human age
        
    elseif i == 2
        % Plot the Vtox=0 control (dotted black)
        plot(age, X(:,1), 'k:', 'LineWidth', 1.6);
    else
        % Plot the fully healthy control (solid green)
        plot(age, X(:,1), '-', 'Color', styleColor, 'LineWidth', 1.8);
    end
    
    % Output the final surviving neuronal fraction at age 100 to the console
    fprintf('%-24s s(100) = %.3f\n', runs{i,1}, interp1(t, X(:,1), 70));
end

%% ========================================================================
% 5. FIGURE FORMATTING
% ========================================================================
% Plot the theoretical resistant neuron floor (s_r = 0.20)
yline(sr, '-.', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2);

% Plot the vertical line indicating the exact age of the catastrophic switch
xline(ageSw, 'r--', 'LineWidth', 1.2);

xlim([30 100]); ylim([0 1.05]);
xlabel('age (yr)');
ylabel('surviving neurons s');
title('Neuronal survival under the age drift');
legend({'full model', 'V_{tox} = 0', 'toxic switch removed', ...
    's_r = 0.20', 'switch'}, 'Location', 'southwest');

% Print the biological switch age to the console
fprintf('microglial switch (full model) at age %.1f\n', ageSw);

%% ========================================================================
% 6. HELPER FUNCTION: LIFETIME ODE SYSTEM
% ========================================================================
% This function evaluates the full 4-ODE system for any given time step.
function dx = lifetime(t, x, sr, d0, da, dtau, Vtox, Ktox, ...
        rm, beta, n, K, gamma, pstim, Amax, cclear, Vc, theta, ...
        ktau, tau0, btau, gtau, dt, p0, alpha)
    
    % Unpack state variables
    s = x(1); m = x(2); a = x(3); tau = x(4);
    
    % Helper equations
    sig = a/(1+a);                  % Saturating microglial uptake
    H = a^n/(K^n + a^n);            % Amyloid-triggered switch probability
    
    % Apply the linear age drift to baseline amyloid production
    p = p0 + alpha*t;
    
    % Calculate the total per-capita death hazard (Phi)
    Phi = d0 + da*sig + dtau*tau^2/(1+tau^2) + Vtox*m^2/(Ktox^2 + m^2);
    
    % 4-ODE System Vector [ds/dt; dm/dt; da/dt; dtau/dt]
    dx = [ -(s-sr)*Phi;
           rm*((1-m)*(beta*H + gamma*m) - m);
           (p + pstim*m)*(1 - a/Amax) - cclear*a - Vc*((1-m)+theta*m)*sig;
           ktau*s*(tau0 + btau*sig + gtau*m) - dt*tau ];
end