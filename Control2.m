clear; clc; close all;

% Neuronal death hazard parameters
sr = 0.20; d0 = 2e-4; da = 0.002; dtau = 0.05;
Vtox0 = 0.12; Ktox = 0.4;  % Baseline microglial toxicity parameters

% Microglial & Switch parameters
rm = 1.0; beta = 4.0; n = 4; sigmaStar = 0.75;
K = sigmaStar/(1-sigmaStar);          % Half-saturation constant (K = 3)
gamma = 0.5;

% Amyloid-beta parameters
pstim0 = 2.5; Amax = 20; cclear = 1.0; Vc = 10; theta0 = 0.3;

% Tau dynamics parameters
ktau = 0.1; tau0 = 0.05; btau = 0.5; gtau = 2.0; dt = 0.1;

% Hill function for the microglial switch
H = @(a) (a.^n)./(K^n + a.^n);

% Analytical m-nullcline (stable resting state of microglia)
mNull = @(a) (-(beta*H(a)+1-gamma) + ...
    sqrt((beta*H(a)+1-gamma).^2 + 4*gamma*beta*H(a)))/(2*gamma);

% Saturating amyloid effect function
sig = @(a) a./(1+a);

% Define a struct array to loop through both the full model and the control
cases = struct( ...
    'name',  {'full model', 'toxic switch removed'}, ...
    'theta', {theta0, 1}, ... % Control: theta=1 means toxic microglia clear just as well as healthy ones
    'pstim', {pstim0, 0}, ... % Control: pstim=0 means toxic microglia do NOT stimulate extra amyloid
    'Vtox',  {Vtox0,  0});    % Control: Vtox=0 means toxic microglia do NOT directly kill neurons

pScan = linspace(2, 60, 581); % Sweep A-beta production

% Setup Figure (Logarithmic Y-axis to visualize the massive hazard jump)
figure('Color','w','Position',[80 80 720 420]);
hold on; box on; set(gca,'YScale','log');

% Loop twice: once for the full model (ic=1), once for the control (ic=2)
for ic = 1:2
    % Extract the specific parameters for the current case
    th = cases(ic).theta;
    ps = cases(ic).pstim;
    Vt = cases(ic).Vtox;
    
    % Define the specific da/dt root-finding function for this case
    F = @(a,p) (p + ps*mNull(a)).*(1-a/Amax) - cclear*a ...
             - Vc*((1-mNull(a))+th*mNull(a)).*sig(a);
         
    nEq = zeros(size(pScan));
    aStar = nan(size(pScan));
    mStar = nan(size(pScan));
    Phi = nan(size(pScan)); % Array to store the equilibrium neuronal hazard
    
    % Loop through every production rate to find equilibria
    for i = 1:numel(pScan)
        p = pScan(i);
        
        % Use a fine grid to locate sign changes for root finding
        agrid = logspace(-2, log10(Amax-0.01), 6000);
        Fa = arrayfun(@(a) F(a,p), agrid);
        roots = [];
        for j = 1:numel(agrid)-1
            if Fa(j)*Fa(j+1) < 0
                roots(end+1) = fzero(@(a) F(a,p), [agrid(j) agrid(j+1)]); 
            end
        end
        nEq(i) = numel(roots);
        
        % If no roots exist (mathematically unlikely here), skip to next p
        if isempty(roots), continue; end
        
        % Always track the lowest stable equilibrium (roots(1)). 
        % When the fold occurs, this root vanishes, forcing the system to jump.
        a = roots(1);
        m = mNull(a);
        
        % Calculate downstream total tau drive (S) based on current (a,m)
        S = tau0 + btau*sig(a) + gtau*m;          
        aStar(i) = a; mStar(i) = m;
        
        % Calculate the total equilibrium neuronal hazard (Phi*)
        Phi(i) = d0 + da*sig(a) + dtau*S^2/(1+S^2) + Vt*m^2/(Ktox^2+m^2);
    end
    
    % Output validation data to the console
    fprintf('%s: max equilibria up to p = 60 is %d\n', cases(ic).name, max(nEq));
    
    if ic == 1
        % Plot the full model hazard curve (Black line, showing the sudden jump)
        plot(pScan, Phi, 'k', 'LineWidth', 1.6);
    else
        % Plot the control hazard curve (Green line, showing smooth progression)
        plot(pScan, Phi, 'Color', [0.15 0.55 0.25], 'LineWidth', 1.8);
        
        % Prove mathematically that the control curve never jumps (monotonic)
        dP = diff(Phi);
        fprintf('Control 2 Phi monotonic: %d\n', all(dP(~isnan(dP)) >= -1e-10));
        fprintf('Control 2 Phi at p = 4, 6.1, 8: %.5f, %.5f, %.5f\n', ...
            interp1(pScan, Phi, 4), interp1(pScan, Phi, 6.1), interp1(pScan, Phi, 8));
    end
end

% Mark the exact mathematical boundaries of the hysteresis window
xline(5.929, 'k:', 'HandleVisibility','off');
xline(6.210, 'k:', 'HandleVisibility','off');

xlim([2 12]); ylim([1e-3 2e-1]);
xlabel('A\beta production p (yr^{-1})');
ylabel('neuronal hazard \Phi^* (yr^{-1})');
title('Control 2: removing the toxic switch abolishes the jump');
legend({'full model','toxic switch removed (\theta=1, p_{stim}=0, V_{tox}=0)'}, ...
    'Location','best');