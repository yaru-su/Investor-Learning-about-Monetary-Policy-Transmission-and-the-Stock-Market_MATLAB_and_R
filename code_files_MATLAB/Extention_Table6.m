%% ========================================================================
%% Table 6 Extension: Asymmetric "J-Shape" Regression Results
%% ========================================================================

% -------------------------------------------------------------------------
% Construct Asymmetric Variables
% -------------------------------------------------------------------------
% Create dummy variables for Tightening (phi > 0) and Easing (phi < 0)
% phi_sq_tight = phi^2 only when phi > 0, else 0
phi_sq_tight = (TT.phi.^2) .* (TT.phi > 0);

% phi_sq_loose = phi^2 only when phi < 0, else 0
phi_sq_loose = (TT.phi.^2) .* (TT.phi < 0);

% -------------------------------------------------------------------------
% Panel A: Market Risk Premium (Data Regression)
% -------------------------------------------------------------------------
% Define Dependent Variable (Data)
RPemp = 12 * TT.LBR_30; 

% Run Regression: RP ~ a + phi_sq_tight + phi_sq_loose
% Note: X1 = a, X2 = Tightening^2, X3 = Easing^2
reg_ext_rp = fitlm([TT.a phi_sq_tight phi_sq_loose], RPemp);

% Compute HAC standard errors
tab_ext_rp = hac_robust(reg_ext_rp);

% Display Results
disp('>>> EXTENSION RESULTS: Risk Premium (Panel A, Data) <<<');
disp('Variables: Intercept | a | Phi^2_Tight | Phi^2_Loose');
disp(tab_ext_rp);

% -------------------------------------------------------------------------
% Panel B: Market Return Volatility (Data Regression)
% -------------------------------------------------------------------------
% Define Dependent Variable (Data)
Volemp = TT.VIX;

% Run Regression: Vol ~ a + phi_sq_tight + phi_sq_loose
reg_ext_vol = fitlm([TT.a phi_sq_tight phi_sq_loose], Volemp);

% Compute HAC standard errors
tab_ext_vol = hac_robust(reg_ext_vol);

% Display Results
disp('>>> EXTENSION RESULTS: Volatility (Panel B, Data) <<<');
disp('Variables: Intercept | a | Phi^2_Tight | Phi^2_Loose');
disp(tab_ext_vol);