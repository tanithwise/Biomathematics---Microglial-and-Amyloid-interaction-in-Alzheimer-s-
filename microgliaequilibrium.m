clear; clc; close all;

%% ========================================================================
% 1. PARAMETER INITIALIZATION (Based on Table 3)
% ========================================================================

% Microglial dynamics parameters
rm = 1.0;           % Microglial transition rate (yr^-1)
beta = 4.0;         % Max amyloid-triggered switch rate
n = 4;              % Hill coefficient for microglial switch (highly sensitive)
sigmaStar = 0.75;   % Saturation threshold for clearance 
K = sigmaStar/(1-sigmaStar);   % Half-saturation constant for switch (K = 3)
gamma = 0.5;        % Cytokine self-amplification rate

% Amyloid-beta dynamics parameters
pstim = 2.5;        % Inflammation-driven A-beta production rate
Amax = 20;          % Carrying capacity / deposition plateau for A-beta
cclear = 1.0;       % Non-microglial clearance rate (yr^-1)
Vc = 10;            % Maximum microglial clearance capacity
theta = 0.3;        % Clearance efficiency penalty for toxic microglia

% Neuronal and Tau dynamics parameters (Downstream)
sr = 0.20;          % Resistant surviving neuron fraction (floor)
d0 = 2e-4;          % Baseline age-related death hazard
da = 0.002;         % A-beta oligomer toxicity coefficient
dtau = 0.05;        % Tangle toxicity coefficient
Vtox = 0.12;        % Maximum microglial neurotoxicity
Ktox = 0.4;         % Half-max toxic fraction for neuronal death
tau0 = 0.05;        % Baseline tau production
btau = 0.5;         % Amyloid-driven tau production
gtau = 2.0;         % Inflammation-driven tau production
dt = 0.1;           % Slow tangle clearance rate (yr^-1)

% Bifurcation Parameter
p = 5.9;  % Baseline A-beta production rate (Swap to 6.1 or 6.3 for other states)

%% ========================================================================
% 2. MATHEMATICAL FUNCTIONS & NULLCLINES
% ========================================================================

% H(a): Hill function governing the amyloid-triggered microglial switch
H = @(a) (a.^n)./(K^n + a.^n);

% mStar(a): The analytical m-nullcline (derived via quadratic formula)
% Solves dm/dt = 0 for m, returning the stable resting state of microglia
mStar = @(a) (-(beta*H(a)+1-gamma) + sqrt((beta*H(a)+1-gamma).^2 + 4*gamma*beta*H(a)))/(2*gamma);

% sig(a): Saturating microglial uptake function (a / 1+a)
sig = @(a) a./(1+a);

% F(a): 1D root-finding function for the whole system. 
% Substitutes the m-nullcline (mStar) into the a-equation to find where da/dt = 0
F = @(a) (p + pstim*mStar(a)).*(1-a/Amax) - cclear*a ...
       - Vc*((1-mStar(a))+theta*mStar(a)).*sig(a);

% Full 2D differential equations (used for the vector field quiver plot)
adot = @(a,m) (p + pstim*m).*(1-a/Amax) - cclear*a - Vc*((1-m)+theta*m).*sig(a);
mdot = @(a,m) rm*((1-m).*(beta*H(a)+gamma*m) - m);

%% ========================================================================
% 3. ROOT FINDING (EQUILIBRIA DETECTION)
% ========================================================================

% Create a fine, logarithmically spaced grid to scan for sign changes.
% This ensures we don't miss closely spaced roots in the bistable window.
agrid = logspace(-2, log10(Amax-0.01), 6000);
Fa = arrayfun(F, agrid);

aEq = []; % Array to store equilibrium A-beta values
for i = 1:length(agrid)-1
    if Fa(i)==0
        aEq(end+1) = agrid(i); % Exact root found on grid
    elseif Fa(i)*Fa(i+1) < 0
        % Sign change detected; use fzero to find the exact root between grid points
        aEq(end+1) = fzero(F,[agrid(i) agrid(i+1)]); 
    end
end

% Calculate the corresponding equilibrium m values
mEq = mStar(aEq);
fprintf('p = %.3f: found %d equilibria\n', p, length(aEq));

%% ========================================================================
% 4. PHASE PLANE PLOTTING (VECTOR FIELD & NULLCLINES)
% ========================================================================

figure; hold on; box on;

% Generate and plot the vector field (quiver plot)
[AA,MM] = meshgrid(0:0.25:5, 0:0.05:1);
dA = adot(AA,MM); 
dM = mdot(AA,MM);
L = sqrt(dA.^2+dM.^2)+eps; % Normalize vector lengths for clean visualization
quiver(AA,MM,dA./L,dM./L,0.5,'Color',[0.7 0.7 0.7], 'HandleVisibility','off');

% Plot the continuous m-nullcline (Orange)
aplot = linspace(0,5,400);
h1 = plot(aplot, mStar(aplot), 'Color',[0.85 0.45 0], 'LineWidth',2);

% Plot the continuous a-nullcline (Blue)
% Derived by solving da/dt = 0 for m analytically
mA = (cclear*aplot + Vc*sig(aplot) - p*(1-aplot/Amax)) ...
   ./ (pstim*(1-aplot/Amax) + Vc*(1-theta)*sig(aplot));
mA(mA<0 | mA>1) = NaN; % Remove mathematically invalid regions outside [0,1]
h2 = plot(aplot, mA, 'b','LineWidth',2);

% Vertical dotted line showing the microglial switch threshold K
xline(K,'r:','LineWidth',1.5, 'HandleVisibility','off');

% Legend setup variables
hStable = [];
hSaddle = [];
legH = [h1 h2];
legN = {'microglia m-nullcline (Eq.3)', 'amyloid a-nullcline (Eq.4)'};

%% ========================================================================
% 5. STABILITY ANALYSIS (JACOBIAN) & HAZARD CALCULATION
% ========================================================================

% Loop through every equilibrium point found by fzero
for k = 1:length(aEq)
    a = aEq(k); 
    m = mEq(k);
    
    % Evaluate Hill function and its derivative at equilibrium
    Hv = H(a); 
    Hp = n*K^n*a^(n-1)/(K^n+a^n)^2;
    
    % Analytically construct the 2x2 (a,m) Jacobian Matrix [J11 J12; J21 J22]
    J11 = -(p+pstim*m)/Amax - cclear - Vc*((1-m)+theta*m)/(1+a)^2;
    J12 = pstim*(1-a/Amax) + Vc*(1-theta)*sig(a);
    J21 = rm*(1-m)*beta*Hp;
    J22 = rm*(-beta*Hv + gamma*(1-2*m) - 1);
    
    Jam = [J11 J12; J21 J22];
    
    % Determine stability: Since trace is strictly negative in this system, 
    % a positive determinant guarantees a stable node, and negative guarantees a saddle.
    isStable = det(Jam) > 0;
    
    % Calculate downstream variables (Tau burden and Neuronal Hazard Phi)
    tau1 = tau0 + btau*sig(a) + gtau*m;
    Phi1 = d0 + da*sig(a) + dtau*tau1^2/(1+tau1^2) + Vtox*m^2/(Ktox^2+m^2);
    taus = sr*tau1;
    Phis = d0 + da*sig(a) + dtau*taus^2/(1+taus^2) + Vtox*m^2/(Ktox^2+m^2);
    
    % Full 4D system eigenvalues (incorporating the s and tau eigenvalues)
    eFull = sort([eig(Jam); -Phis; -dt]);
    
    typ = "stable node";
    if ~isStable, typ = "saddle"; end
    
    % Print analytical results to the console
    fprintf('eq%d: a*=%.3f m*=%.3f Phi(s=1)=%.4f Phi(s=sr)=%.4f %s\n', k, a, m, Phi1, Phis, typ);
    
    % Plot the equilibrium point markers dynamically
    if isStable
        % Stable nodes are plotted as solid black dots
        h = plot(a, m, 'ko', 'MarkerFaceColor','k', 'MarkerSize',9);
        if isempty(hStable)
            hStable = h;
            legH(end+1) = h; 
            legN{end+1} = 'stable node'; 
        else
            set(h, 'HandleVisibility','off'); % Prevent duplicate legend entries
        end
    else
        % Unstable saddle nodes are plotted as hollow white dots
        h = plot(a, m, 'ko', 'MarkerFaceColor','w', 'MarkerSize',9, 'LineWidth',1.5);
        if isempty(hSaddle)
            hSaddle = h;
            legH(end+1) = h; 
            legN{end+1} = 'saddle'; 
        else
            set(h, 'HandleVisibility','off'); % Prevent duplicate legend entries
        end
    end
end

%% ========================================================================
% 6. FIGURE FORMATTING
% ========================================================================

xlim([0 5]); ylim([0 1]);
xlabel('A\beta load a'); 
ylabel('toxic microglia m');
title(sprintf('a-m phase plane at p = %.1f', p));
legend(legH, legN, 'Location','best');