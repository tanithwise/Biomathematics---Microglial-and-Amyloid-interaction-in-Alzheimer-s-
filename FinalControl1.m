% control1_linear_clearance.m
clear; clc; close all;

% Microglial & Switch Parameters
rm = 1.0; beta = 4.0; n = 4; sigmaStar = 0.75;
K = sigmaStar/(1-sigmaStar);          % Half-saturation constant (K = 3)
gamma = 0.5;

% Amyloid-beta Parameters
pstim = 2.5; Amax = 20; cclear = 1.0; Vc = 10; theta = 0.3;

% H(a): Hill function governing the highly sensitive microglial switch
H = @(a) (a.^n)./(K^n + a.^n);

% mNull(a): Analytical m-nullcline (stable resting state of microglia)
mNull = @(a) (-(beta*H(a)+1-gamma) + ...
    sqrt((beta*H(a)+1-gamma).^2 + 4*gamma*beta*H(a)))/(2*gamma);

% Fof: Wrapper function to evaluate da/dt = 0. 
% The 'mode' argument acts as a switch between 'sat' (full model) and 'lin' (control).
% eff = 1-(1-theta)*m = (1-m)+theta*m
Fof = @(a,p,mode) localF(a,p,mode,mNull,pstim,Amax,cclear,Vc,theta);

pScan = linspace(2, 60, 581); % Sweep A-beta production up to 60 (10x past the fold)

% Preallocate arrays to store the number of roots and their values
nSat = zeros(size(pScan));
nLin = zeros(size(pScan));
aSat = nan(numel(pScan), 3);  % Saturating model can have up to 3 roots (bistability)
aLin = nan(numel(pScan), 1);  % Linear model should only ever have 1 root

% Loop through every production rate and find the equilibrium points
for i = 1:numel(pScan)
    p = pScan(i);
    
    % Find equilibrium points for both models
    rS = rootsOf(@(a) Fof(a,p,'sat'));
    rL = rootsOf(@(a) Fof(a,p,'lin'));
    
    nSat(i) = numel(rS);
    nLin(i) = numel(rL);
    
    aSat(i,1:numel(rS)) = rS;
    if ~isempty(rL)
        aLin(i) = rL(1);      % Store the single root for the linear model
    end
end

fprintf('Saturating clearance: max equilibria in scan = %d\n', max(nSat));
fprintf('Linear clearance:     max equilibria up to p = 60 is %d\n', max(nLin));

% Prove that the saddle-node fold never occurs in the control model
if max(nLin) <= 1
    fprintf('No saddle-node fold for any p up to 60.\n');
else
    fprintf('Fold detected under linear clearance — check the grid.\n');
end

% Verify the linear control amyloid load rises smoothly without jumping
dA = diff(aLin);
fprintf('Linear a* monotonic nondecreasing: %d\n', all(dA >= -1e-8));
fprintf('Linear a* at p = 4, 6.1, 8, 60: %.3f, %.3f, %.3f, %.3f\n', ...
    interp1(pScan,aLin,4), interp1(pScan,aLin,6.1), ...
    interp1(pScan,aLin,8), aLin(end));

% Calculate the Determinant of the Jacobian for the linear model across the scan.
% A saddle-node bifurcation requires det(J) = 0. This loop proves it never hits zero.
detLin = nan(size(pScan));
for i = 1:numel(pScan)
    if isnan(aLin(i)), continue; end
    detLin(i) = detJam(aLin(i), mNull(aLin(i)), pScan(i), 'lin', ...
        rm, beta, n, K, gamma, pstim, Amax, cclear, Vc, theta, H);
end
fprintf('Linear det Jam range on the scan: [%.3f, %.3f] (never zero)\n', ...
    min(detLin), max(detLin));


figure('Color','w','Position',[80 80 720 420]);
hold on; box on;

% Separate the 3 roots of the saturating model into stable nodes and the saddle
stableSat = aSat;
saddleSat = nan(size(pScan));
for i = 1:numel(pScan)
    if nSat(i) == 3
        saddleSat(i) = aSat(i,2); % The middle root is the unstable saddle
        stableSat(i,2) = nan;     % Remove the saddle from the stable array
    end
end

% Plot the full model (black stable nodes, red unstable saddle forming the S-curve)
plot(pScan, stableSat(:,1), 'k', 'LineWidth', 1.6);
plot(pScan, stableSat(:,3), 'k', 'LineWidth', 1.6, 'HandleVisibility','off');
plot(pScan, saddleSat, 'r', 'LineWidth', 1.4);

% Plot the linear clearance control (smooth, monotonically increasing grey line)
plot(pScan, aLin, 'Color', [0.45 0.45 0.45], 'LineWidth', 1.8);

% Mark the hysteresis window boundaries (p1 and p2)
xline(5.929, 'k:', 'HandleVisibility','off');
xline(6.210, 'k:', 'HandleVisibility','off');

% Formatting
xlim([2 20]); ylim([0 8]);          % Window is zoomed to show the fold
xlabel('A\beta production p (yr^{-1})');
ylabel('equilibrium amyloid a^*');
title('Control 1: linear clearance removes the fold');
legend({'saturating clearance','saddle (saturating)','linear clearance'}, ...
    'Location','northwest');


% localF: Evaluates da/dt = 0, toggling between saturating and linear microglial uptake
function F = localF(a, p, mode, mNull, pstim, Amax, cclear, Vc, theta)
    m = mNull(a);
    eff = (1-m) + theta*m;
    
    if strcmp(mode,'sat')
        uptake = a./(1+a); % Full model uses saturable Hill function
    else
        uptake = a;        % Control 1 uses strictly linear uptake
    end
    
    F = (p + pstim*m).*(1 - a/Amax) - cclear*a - Vc*eff.*uptake;
end

% -------------------------------------------------------------------------
% detJam: Calculates the determinant of the (a,m) Jacobian matrix.
% J11 and J12 are mathematically different depending on the clearance mode.
function d = detJam(a, m, p, mode, rm, beta, n, K, gamma, pstim, Amax, cclear, Vc, theta, H)
    eff = (1-m) + theta*m;
    
    if strcmp(mode,'sat')
        J11 = -(p+pstim*m)/Amax - cclear - Vc*eff/(1+a)^2;
        J12 = pstim*(1-a/Amax) + Vc*(1-theta)*a/(1+a);
    else
        J11 = -(p+pstim*m)/Amax - cclear - Vc*eff;          % Derivative of linear 'a' is a constant
        J12 = pstim*(1-a/Amax) + Vc*(1-theta)*a;
    end
    
    Hp = n*K^n*a^(n-1)/(K^n + a^n)^2;
    J21 = rm*(1-m)*beta*Hp;
    J22 = rm*(-beta*H(a) + gamma*(1-2*m) - 1);
    
    d = J11*J22 - J12*J21; % Determinant calculation
end

% -------------------------------------------------------------------------
% rootsOf: Locates roots by scanning a fine logarithmic grid for sign changes, 
% then using fzero to pinpoint the exact equilibrium coordinate.
function roots = rootsOf(F)
    agrid = logspace(-2, log10(19.5), 8000);
    Fa = arrayfun(F, agrid);
    roots = [];
    for i = 1:numel(agrid)-1
        if Fa(i)*Fa(i+1) < 0
            roots(end+1) = fzero(F, [agrid(i) agrid(i+1)]); 
        end
    end
    roots = unique(round(roots, 6));   % Drop duplicates from a flat crossing
end