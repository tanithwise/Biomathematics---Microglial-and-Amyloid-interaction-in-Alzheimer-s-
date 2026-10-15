clear; clc;

% These parameters match the 'Complete parameter set' (Table 3).
rm = 1.0; beta = 4.0; n = 4; sigmaStar = 0.75;
K = sigmaStar/(1-sigmaStar); % Half-saturation constant (K = 3)
gamma = 0.5;
pstim = 2.5; Amax = 20; cclear = 1.0; Vc = 10; theta = 0.3;

% Fixed production rate set to the bistable window to analyze the healthy node
p = 6.1;

H  = @(a) (a.^n)./(K^n + a.^n);                             % Hill function for switch
Hp = @(a) n*K^n*a.^(n-1)./((K^n + a.^n).^2);                % First derivative of the Hill function, H'(a)
sig = @(a) a./(1+a);                                        % Saturating microglial uptake

% The 2D subsystem differential equations
adot = @(a,m) (p + pstim*m).*(1-a/Amax) - cclear*a - Vc*((1-m)+theta*m).*sig(a);
mdot = @(a,m) rm*((1-m).*(beta*H(a)+gamma*m) - m);

% Using the previously calculated healthy equilibrium coordinates at p = 6.1
a = 1.1927; m = 0.1455;

% Calculate the four elements of the Jacobian matrix analytically based on 
% the hand-derived partial derivatives (Eqs. 5-8 in the methodology).
J11 = -(p+pstim*m)/Amax - cclear - Vc*((1-m)+theta*m)/(1+a)^2;      % da_dot/da
J12 = pstim*(1-a/Amax) + Vc*(1-theta)*sig(a);                       % da_dot/dm
J21 = rm*(1-m)*beta*Hp(a);                                          % dm_dot/da
J22 = rm*(-beta*H(a) + gamma*(1-2*m) - 1);                          % dm_dot/dm

Jam_analytic = [J11 J12; J21 J22]; % Construct the 2x2 analytical matrix

% This section numerically approximates the Jacobian to prove that the 
% hand-derived analytical equations above contain no mathematical errors.

rhs = @(x) [adot(x(1),x(2)); mdot(x(1),x(2))]; % Vectorize the ODEs
x0 = [a; m];                                   % Point of evaluation
h = 1e-6;                                      % Extremely small step size for precision

Jnum = zeros(2,2); % Initialize the numerical matrix

% Loop through variables a and m to calculate partial derivatives
for j = 1:2
    dx = [0;0]; 
    dx(j) = h;
    % Central difference formula: [f(x + h) - f(x - h)] / 2h
    Jnum(:,j) = (rhs(x0+dx) - rhs(x0-dx))/(2*h);
end

tr = trace(Jam_analytic);   % Trace = J11 + J22
dtm = det(Jam_analytic);    % Determinant = (J11*J22) - (J12*J21)

% Output the exact matrices to command window
fprintf('Jam analytic =\n'); disp(Jam_analytic);
fprintf('Jam numeric =\n'); disp(Jnum);

% Output stability criteria. 
% For a stable node, trace must be < 0 and determinant must be > 0.
fprintf('tr = %.4f, det = %.4f\n', tr, dtm);

% Compare the two matrices. A very small error (e.g., 10^-10) proves the 
% analytical partial derivatives are perfectly accurate.
fprintf('max |analytic-numeric| = %.2e\n', max(abs(Jam_analytic-Jnum),[],'all'));

% Print the actual eigenvalues to confirm stability analytically
fprintf('eig = %.4f, %.4f\n', eig(Jam_analytic));