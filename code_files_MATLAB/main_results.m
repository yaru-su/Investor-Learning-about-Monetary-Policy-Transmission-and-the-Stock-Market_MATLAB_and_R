%% 2_main_results

%% 1. Figure 1: GDP Growth, Output Gap, Inflation, and Federal Funds Rate

figure(1);
% First panel: Real GDP Growth
subplot(4,1,1);
plot(TT.date, TT.GDPGrowth);
ylabel('GDP growth');
% Second panel: Output Gap
subplot(4,1,2);
plot(TT.date, TT.OutputGap);
ylabel('Output gap');
% Third panel: Inflation Rate
subplot(4,1,3);
plot(TT.date, TT.inflation);
ylabel('Inflation');
% Fourth panel: Federal Funds Rate
subplot(4,1,4);
plot(TT.date, TT.rfnom);
ylabel('Fed funds rate');
% Save the figure as a PNG
saveas(gcf, 'Figure1_MacroVariables.png');


%% 2. Consumption Growth: ML Estimation of (mudeltabar,sigmadelta) mudeltabar not needed in model
% Requires: LikelihoodFunc_Delta, GetVarMatrixParam_Delta

% Set optimization options for fminsearch
options = optimset('TolX',1e-4,'TolFun',1e-4,'Algorithm','interior-point');   
% Initial guess for the parameters [mudeltabar; sigmadelta]
paramdelta0=[0.02;0.02];
% Maximize the likelihood function to estimate sigmadelta
[paramdelta,fval]=fminsearch(@(x)LikelihoodFunc_Delta(TT.GDPGrowth,x),paramdelta0,options);
paramdelta;

% Compute variance-covariance matrix of estimated parameters
[Vdelta]=GetVarMatrixParam_Delta(TT.GDPGrowth,paramdelta);
% Compute standard errors from the diagonal of the variance-covariance matrix
paramdelta_std = sqrt(diag(Vdelta));
% Number of observations
nobs = length(TT.GDPGrowth); 
% Compute t-stat of estimated parameters
paramdelta_tstats  = abs(paramdelta./paramdelta_std);
% Compute two-sided p-values
paramdelta_pvalues = 2*(1-tcdf(paramdelta_tstats,nobs-numel(paramdelta)));

% Display estimated parameters and their statistics
T = table(paramdelta, paramdelta_std, paramdelta_tstats, paramdelta_pvalues, ...
    'VariableNames', {'Estimate','Std_Error','tStat','pValue'});
disp(T)


%% 3. Output Gap: ML estimation of (sigmay,lambday)
% Requires: LikelihoodFuncy, GetVarMatrixParamy

% Estimate an initial guess for λ_y using AR(1) regression
reg=nwest(TT.OutputGap(2:end,1),[ones(length(TT.OutputGap)-1,1) TT.OutputGap(1:end-1,1)],0);
% Convert AR(1) coefficient into continuous-time mean-reversion rate
lamby_reg=-log(reg.beta(2))/Delta;

% Set initial parameter vector (σ_y, λ_y)
paramy0=[std(TT.OutputGap);lamby_reg];
% Maximize the likelihood function to estimate parameters
[paramy,fval]=fminsearch(@(x)LikelihoodFuncy(TT.OutputGap,x),paramy0,options);

% Compute variance-covariance matrix of estimated parameters
[Vy]=GetVarMatrixParamy(TT.OutputGap,paramy);
% Compute standard errors from the diagonal of the variance-covariance matrix
paramy_std = sqrt(diag(Vy));
% Number of observations
nobs = length(TT.OutputGap);
% Compute t-stat of estimated parameters
paramy_tstats  = abs(paramy./paramy_std);
% Compute two-sided p-values
paramy_pvalues = 2*(1-tcdf(paramy_tstats,nobs-numel(paramy)));

% Display estimated parameters and their statistics
ParamNames = {'mu_delta_bar','sigma_delta'}';
paramy         = paramy(:);
paramy_std     = paramy_std(:);
paramy_tstats  = paramy_tstats(:);
paramy_pvalues = paramy_pvalues(:);

ResultsTable_OutputGap = table(ParamNames, paramy, paramy_std, ...
    paramy_tstats, paramy_pvalues, ...
    'VariableNames', {'Parameter','Estimate','StdError','tStat','pValue'});
disp(ResultsTable_OutputGap);

%% 4. Taylor Rule: ML Estimation of (rNbar,betapi,betay,volr) volr not needed in model
% Requires: LikelihoodFunc_Rate, GetVarMatrixParam_Rate

% Define initial parameter guesses
paramr0=[mean(TT.rfnom);0.25;0.15;0.02];
% Use fminsearch to perform Maximum Likelihood estimation
[paramr,fval]=fminsearch(@(x)LikelihoodFunc_Rate(TT.rfnom,TT.inflation,TT.OutputGap,x),paramr0,options);
% Compute log-likelihood and model-implied variables
[sumloglik,logLik,rfTaylor,stdeps]=LikelihoodFunc_Rate(TT.rfnom,TT.inflation,TT.OutputGap,paramr);

% Compute variance-covariance matrix of estimated parameters
[Vr]=GetVarMatrixParam_Rate(TT.rfnom,TT.inflation,TT.OutputGap,paramr);
% Compute standard errors from the diagonal of the variance-covariance matrix
paramr_std = sqrt(diag(Vr));
% Number of observations
nobs = length(TT.inflation); 
% Compute t-stat of estimated parameters
paramr_tstats  = abs(paramr./paramr_std);
% Compute two-sided p-values
paramr_pvalues = 2*(1-tcdf(paramr_tstats,nobs-numel(paramr)));

% Extract estimated rNbar
rNbar_estimated=paramr(1);

% Display estimated parameters and their statistics
disp('Estimated Parameters for Taylor Rule:');
disp(paramr);
disp('Standard Errors for Taylor Rule:');
disp(paramr_std);
disp('T-Statistics for Taylor Rule:');
disp(paramr_tstats);
disp('P-Values for Taylor Rule:');
disp(paramr_pvalues);


%% 5. Inflation and Transmission Mechanism: ML estimation of inflation parameters (sigmapi,pibar,lambdapi,sigmaa,lambdaa)
% Requires: LikelihoodFunc, GetVarMatrixParam

% Compute sample covariance between current and lagged inflation
COV=cov(TT.inflation(2:end,1),TT.inflation(1:end-1,1));
% Compute initial mean-reversion speed λπ
lambdapi0=-1/Delta*log(COV(1,2)/COV(1,1));
% Compute initial volatility σπ
sigmapi0=sqrt(var(TT.inflation)*2*lambdapi0);
% Compute initial long-run inflation target π̄
pibar0=mean(TT.inflation);

% Combine all initial guesses into parameter vector
param0=[sigmapi0;pibar0;0.1;0.2;0.5];
% Estimate parameters via Maximum Likelihood
[param,fval]=fminsearch(@(x)LikelihoodFunc(TT.inflation,TT.fedfund,x),param0,options);
% Recompute likelihood and filtered states using estimated parameters
[sumloglik,logLik,a,nua,LTMpi,stdeps]=LikelihoodFunc(TT.inflation,TT.fedfund,param);

% Compute variance-covariance matrix of parameters
[V]=GetVarMatrixParam(TT.inflation,TT.fedfund,param);
% Compute standard errors from the diagonal of the variance-covariance matrix
param_std = sqrt(diag(V));
% Number of observations
nobs = length(TT.inflation); 
% Compute t-stat of estimated parameters
param_tstats  = abs(param./param_std);
% Compute two-sided p-values
param_pvalues = 2*(1-tcdf(param_tstats,nobs-numel(param)));

% Display estimated parameters and their statistics for inflation
disp('Estimated Parameters for Inflation:');
disp(param);
disp('Standard Errors for Inflation:');
disp(param_std);
disp('T-Statistics for Inflation:');
disp(param_tstats);
disp('P-Values for Inflation:');
disp(param_pvalues);

% Insert ML-Kalman estimated state variables in timetable
TT.a=a;  % learning (transmission) coefficient
TT.nua=nua;  % uncertainty about a_t
TT.LTMpi=LTMpi;  % long-term mean of inflation
TT.phi=rfTaylor-rNbar_estimated;  % monetary policy deviation

% Analyze properties of estimated state variables
% (a) Relationship between changes in a_t and inflation
change_a=diff(TT.a);
change_pi=diff(TT.inflation);
dummy=(TT.inflation>mean(TT.inflation));
corr(change_a.*dummy(2:end),change_pi.*dummy(2:end));
corr(change_a.*(1-dummy(2:end)),change_pi.*(1-dummy(2:end)));

% (b) Moments of φ = r_N,t - r̄_N (monetary policy deviation)
regautocorrphi=nwest(TT.phi(2:end),[ones(length(TT.phi(2:end)),1) TT.phi(1:end-1)],0);
autocorrphi=regautocorrphi.beta(2);
Volphi=std(TT.phi);

% (c) Moments of long-term mean of inflation (LTMπ)
reg=nwest(LTMpi(2:end,1),[ones(length(LTMpi(2:end,1)),1) LTMpi(1:end-1,1)],0);
autocorLTMpi=reg.beta(2);
volLTMpi=std(LTMpi);

% (d) Correlation with 5-year inflation expectations (footnote 11)
corr(TT.CPI5YR,TT.LTMpi,"rows","complete")
reg=fitlm(TT.LTMpi,TT.CPI5YR);
hac_robust(reg)


%% 6. Figure 2: Investors' Estimate of the Transmission Coefficient

figure(2);
% Left y-axis: transmission coefficient
yyaxis left
plot(TT.date, TT.a, 'b'); % blue line
ylabel('Transmission coefficient $\widehat{a}_t$', 'Interpreter', 'latex')

% Right y-axis: inflation
yyaxis right
plot(TT.date, TT.inflation, 'r'); % red line
ylabel('Inflation $\pi_t$', 'Interpreter', 'latex')


%% 7. Figure 3: Monetary Policy Stance and Long-Term Inflation Drift

figure(3);
% First panel
subplot(2,1,1);
plot(TT.date, TT.phi);
ylabel('Monetary policy stance $\phi_t$', 'Interpreter', 'latex');
% Second panel
subplot(2,1,2);
plot(TT.date, TT.LTMpi);
ylabel('Long-term inflation drift $\widehat{\pi}_t^A$', 'Interpreter', 'latex');


%% 8. Table 3: Parameter Values Estimated by Maximum Likelihood

% sigmadelta,sigmay,lambday,rNbar,betapi,betay,sigmapi,pibar,lambdapi,sigmaa,lambdaa
ParamEstimTable=[paramdelta(2,1) paramdelta_std(2,1);
    paramy paramy_std;
    paramr(1:end-1,1) paramr_std(1:end-1,1);
    param param_std]


%% 9. SAVE ML-KALMAN ESTIMATED STATE VARIABLES

% saves the state variables estimated from the ML–Kalman filter into a CSV file for table 7-9
writetimetable(TT(:,{'GDPGrowth','rfnom','inflation','a','nua','phi','OutputGap','LTMpi'}),'StateVarsOct2024.csv','Delimiter',',')


%% 10. S&P500 real yearly log dividend growth (used to calibrate sigmaD in model)

divgrowth=log(TT.d12(13:end)./TT.d12(1:end-12))-TT.inflation(13:end);
std(divgrowth);


%% 11. Table 4: Asset-pricing Moments

%(a) Data
%real and nominal risk-free rate
meanrfreal=mean(TT.rfreal)
meanrfnom=mean(TT.rfnom)

%S&P500 excess returns
meanSP=12*mean(TT.exret)
volSP=sqrt(12)*std(TT.exret)


%(b) Table 4, Model is obtained in Mathematica by setting pi=pibar and a=phi=0
%Model-implied time series
DataModel=readtable('data/ModelImpliedTimeSeriesNov2024'); %these series come from Mathematica
DataModel.date=TT.date;
TTmodel=table2timetable(DataModel);


%% 12. Table 5: Expected Output Growth, Real Interest Rate, and Log Price-Dividend Ratio (X) vs. Inflation and the Output Gap (Y).

%% (a) Col 1,2: Expected output growt

% -------------------------------------------------------------------------
% Model regression
% -------------------------------------------------------------------------
% Run OLS regression: regress model-implied expected output growth on inflation and output gap.
reg=fitlm([TT.inflation TT.OutputGap],TTmodel.mudelta);
% Compute HAC-robust (Newey–West) standard errors for time series data.
tab_model=hac_robust(reg);
% Identify which observations (rows) were used (exclude missing values).
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
% Compute standardized effect of a one-standard-deviation increase in inflation and output gap on expected output growth.
OneStdIncrease_x=std([TT.inflation(usedObs) TT.OutputGap(usedObs)]).*reg.Coefficients.Estimate(2:3)';
% Compute the standard deviation of the dependent variable (Y).
OneStd_y=std(TTmodel.mudelta(usedObs));
% Express the impact as a ratio: ΔY(σx)/σy
OneStdIncrease_x/OneStd_y;

% -------------------------------------------------------------------------
% Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.inflation TT.OutputGap],TT.mudelta);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
OneStdIncrease_x=std([TT.inflation(usedObs) TT.OutputGap(usedObs)]).*reg.Coefficients.Estimate(2:3)';
OneStd_y=std(TT.mudelta(usedObs));
OneStdIncrease_x/OneStd_y;

%Table 5, col 1,2
tab_model
tab_data


%% (b) Col 3,4: real risk-free rate

% -------------------------------------------------------------------------
% Model regression
% -------------------------------------------------------------------------
reg=fitlm([TT.inflation TT.OutputGap],TTmodel.rfreal);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
OneStdIncrease_x=std([TT.inflation(usedObs) TT.OutputGap(usedObs)]).*reg.Coefficients.Estimate(2:3)';
OneStd_y=std(TTmodel.rfreal(usedObs));
OneStdIncrease_x/OneStd_y;

% -------------------------------------------------------------------------
% Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.inflation TT.OutputGap],TT.rfreal);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
OneStdIncrease_x=std([TT.inflation(usedObs) TT.OutputGap(usedObs)]).*reg.Coefficients.Estimate(2:3)';
OneStd_y=std(TT.rfreal(usedObs));
OneStdIncrease_x/OneStd_y;

% Table 5, col 3,4
tab_model
tab_data


%% (c) Col 5,6: log PD ratio

% -------------------------------------------------------------------------
% Model regression
% -------------------------------------------------------------------------
reg=fitlm([TT.inflation TT.OutputGap],TTmodel.pd);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
OneStdIncrease_x=std([TT.inflation(usedObs) TT.OutputGap(usedObs)]).*reg.Coefficients.Estimate(2:3)';
OneStd_y=std(TTmodel.pd(usedObs));
OneStdIncrease_x/OneStd_y;

% -------------------------------------------------------------------------
% Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.inflation TT.OutputGap],TT.pd);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
OneStdIncrease_x=std([TT.inflation(usedObs) TT.OutputGap(usedObs)]).*reg.Coefficients.Estimate(2:3)';
OneStd_y=std(TT.pd(usedObs));
OneStdIncrease_x/OneStd_y;

%Table 5, col 5,6
tab_model
tab_data


%% 13. Table 6:  Market Risk Premium and Return Volatility (X) vs. State Variables (Y)

%% Panel A: market risk premium

% -------------------------------------------------------------------------
% Panel A, col 1,2: Model regression
% -------------------------------------------------------------------------
% corr between Martin(model) and ChabiYo RPs(data)
corr(TT.RP1, 12*TT.LBR_30, 'Rows', 'complete'); % 'Rows','complete' means ignore rows with missing values (NaN).
% Define risk premium variables
RPmodel=TTmodel.RP; % Model-implied risk premium from theoretical model (TTmodel)
RPemp=12*TT.LBR_30; %annualized risk premium (1-month-ahead, preference-free) ChabiYo

% Fit a linear model (OLS): RPmodel_t = α + β * a_t + ε_t
reg=fitlm([TT.a],RPmodel);
% Compute HAC-robust (Newey–West) standard errors
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

% -------------------------------------------------------------------------
% Panel A, col 1,2: Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a],RPemp);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

%Table 6, Panel A, col 1,2
tab_model
tab_data


% -------------------------------------------------------------------------
% Panel A, col 3,4: Model regression
% -------------------------------------------------------------------------
reg=fitlm([TT.phi.^2],RPmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

% -------------------------------------------------------------------------
% Panel A, col 3,4: Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.phi.^2],RPemp);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

%Table 6, Panel A, col 3,4
tab_model
tab_data


% -------------------------------------------------------------------------
% Panel A, col 5,6: Model regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a TT.phi.^2],RPmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

% -------------------------------------------------------------------------
% Panel A, col 5,6: Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a TT.phi.^2],RPemp);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

%Table 6, Panel A, col 5,6
tab_model
tab_data


% -------------------------------------------------------------------------
% Panel A, col 7,8: Model regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a TT.phi TT.phi.^2],RPmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
OneStdIncrease_x=std([[TT.a(usedObs) TT.phi(usedObs) TT.phi(usedObs).^2]]).*reg.Coefficients.Estimate(2:4)';
OneStd_y=std(RPmodel(usedObs));
OneStdIncrease_x/OneStd_y;

% -------------------------------------------------------------------------
% Panel A, col 7,8: Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a TT.phi TT.phi.^2],RPemp);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
OneStdIncrease_x=std([[TT.a(usedObs) TT.phi(usedObs) TT.phi(usedObs).^2]]).*reg.Coefficients.Estimate(2:4)';
OneStd_y=std(RPemp(usedObs));
OneStdIncrease_x/OneStd_y;

%Table 6, Panel A, col 7,8
tab_model
tab_data


%% Panel B: Volatility

Volmodel=TTmodel.vol;
Volemp=TT.VIX;

% -------------------------------------------------------------------------
% Panel B, col 1,2: Model regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a],Volmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

% -------------------------------------------------------------------------
% Panel B, col 1,2: Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a],Volemp);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

%Table 6, Panel B, col 1,2
tab_model
tab_data


% -------------------------------------------------------------------------
% Panel B, col 3,4: Model regression
% -------------------------------------------------------------------------
reg=fitlm([TT.phi.^2],Volmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

% -------------------------------------------------------------------------
% Panel B, col 3,4: Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.phi.^2],Volemp);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

%Table 6, Panel B, col 3,4
tab_model
tab_data


% -------------------------------------------------------------------------
% Panel B, col 5,6: Model regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a TT.phi.^2],Volmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

% -------------------------------------------------------------------------
% Panel B, col 5,6: Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a TT.phi.^2],Volemp);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

%Table 6, Panel B, col 5,6
tab_model
tab_data


% -------------------------------------------------------------------------
% Panel B, col 7,8: Model regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a TT.phi TT.phi.^2],Volmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
OneStdIncrease_x=std([[TT.a(usedObs) TT.phi(usedObs) TT.phi(usedObs).^2]]).*reg.Coefficients.Estimate(2:4)';
OneStd_y=std(Volmodel(usedObs));
OneStdIncrease_x/OneStd_y;

% -------------------------------------------------------------------------
% Panel B, col 7,8: Data regression
% -------------------------------------------------------------------------
reg=fitlm([TT.a TT.phi TT.phi.^2],Volemp);
tab_data=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
OneStdIncrease_x=std([[TT.a(usedObs) TT.phi(usedObs) TT.phi(usedObs).^2]]).*reg.Coefficients.Estimate(2:4)';
OneStd_y=std(Volemp(usedObs));
OneStdIncrease_x/OneStd_y;

%Table 6, Panel B, col 7,8
tab_model
tab_data


%% 14. Table 7: Market Risk Premium and Return Volatility (X) vs. Model-Free Proxies (Y)

%% Panel A: Market risk premium vs. model-free proxies
% phi_proxy is inflation deviations from its mean; 'omit' ignores NaNs when computing the mean.
phi_proxy=TT.inflation-mean(TT.inflation,'omit');
% Pearson correlation between TT.phi and phi_proxy; use only rows with no NaNs in either series.
corr(TT.phi,phi_proxy,'rows','complete')

%compute long-term inflation expectation: mean across horizons 5 to 30
ExpectedInflationMat = [];
for i = 5:1:30 % i = [5, 6, 7, ..., 30]
    varname = sprintf('x%dYearExpectedInflation', i);
    ExpectedInflationMat = [ExpectedInflationMat, TT.(varname)]; % Horizontally concatenate each series into the matrix (aligned by TT's time index).
end
% Compute the mean of the expected inflation matrix across the specified horizons
LTExpectedInflationMean = mean(ExpectedInflationMat, 2);
% 5th percentile of inflation as a "low-inflation threshold".
LowQuantileInflation=quantile(TT.inflation,0.05);
% Build a_proxy. Denominator uses elementwise max(low-quantile threshold, actual inflation) to avoid explosive ratios at very low inflation.
a_proxy=-LTExpectedInflationMean./max(LowQuantileInflation,TT.inflation); 
% Correlate model-implied a (TT.a) with the proxy a_proxy using complete cases.
corr(TT.a,a_proxy,'rows','complete')

% Table 7, Panel A, col. 1
% OLS regression (with intercept by default): RPemp_t = β0 + β1 * a_proxy_t + ε_t.
reg=fitlm(a_proxy,RPemp);
% Compute HAC-robust (Newey–West) standard errors
tab_data=hac_robust(reg)
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

% Table 7, Panel A, col. 2
reg=fitlm(phi_proxy.^2,RPemp);
tab_data=hac_robust(reg)

% Table 7, Panel A, col. 3
reg=fitlm([a_proxy phi_proxy.^2],RPemp);
tab_data=hac_robust(reg)

% Table 7, Panel A, col. 4
reg=fitlm([a_proxy phi_proxy phi_proxy.^2],RPemp);
tab_data=hac_robust(reg)


%% Panel B: Market return volatility vs. model-free proxies
% Table 7, Panel B, col. 1
reg=fitlm(a_proxy,Volemp);
tab_data=hac_robust(reg)
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
usedObsBench=usedObs;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

% Table 7, Panel B, col. 2
reg=fitlm(phi_proxy(usedObsBench).^2,Volemp(usedObsBench));
tab_data=hac_robust(reg)

% Table 7, Panel B, col. 3
reg=fitlm([a_proxy phi_proxy.^2],Volemp);
tab_data=hac_robust(reg)

% Table 7, Panel B, col. 4
reg=fitlm([a_proxy phi_proxy phi_proxy.^2],Volemp);
tab_data=hac_robust(reg)


%% 15. Table 8: Impact of Monetary Uncertainty on Market Risk Premium and Volatility

nua_proxy=TT.MPU1; %monetary policy uncertainty, Bloom
IndTime=~isnan(nua_proxy); % IndTime marks the non-missing (valid) observations.

% Define time periods with high monetary policy uncertainty.
Times_HighnuaTrue=(TT.nua>quantile(TT.nua(~isnan(TT.nua)),0.5));
Times_Highnua=double(nua_proxy>quantile(nua_proxy(~isnan(nua_proxy)),0.5));
Times_Highnua(isnan(nua_proxy)) = NaN; % Set NaN where nua_proxy is NaN

%% Panel A: State Variables

% Table 8, Panel A, col. 1
reg=fitlm([TT.a TT.phi TT.phi.^2 Times_HighnuaTrue],RPemp);
tab_data=hac_robust(reg)

% Table 8, Panel A, col. 2
reg=fitlm([TT.a TT.phi TT.phi.^2 Times_HighnuaTrue (TT.phi.^2).*Times_HighnuaTrue],RPemp);
tab_data=hac_robust(reg)

% Table 8, Panel A, col. 3
reg=fitlm([TT.a TT.phi TT.phi.^2 Times_HighnuaTrue],Volemp);
tab_data=hac_robust(reg)

% Table 8, Panel A, col. 4
reg=fitlm([TT.a TT.phi TT.phi.^2 Times_HighnuaTrue (TT.phi.^2).*Times_HighnuaTrue],Volemp);
tab_data=hac_robust(reg)


%% Panel B: Model-Free Proxies

% Table 8, Panel B, col. 1
reg=fitlm([a_proxy phi_proxy phi_proxy.^2 Times_Highnua],RPemp);
tab_data=hac_robust(reg)

% Table 8, Panel B, col. 2
reg=fitlm([a_proxy phi_proxy phi_proxy.^2 Times_Highnua (phi_proxy.^2).*Times_Highnua],RPemp);
tab_data=hac_robust(reg)

% Table 8, Panel B, col. 3
reg=fitlm([a_proxy phi_proxy phi_proxy.^2 Times_Highnua],Volemp);
tab_data=hac_robust(reg)

% Table 8, Panel B, col. 4
reg=fitlm([a_proxy phi_proxy phi_proxy.^2 Times_Highnua (phi_proxy.^2).*Times_Highnua],Volemp);
tab_data=hac_robust(reg)


%% 16. Predictive regressions for future realized excess returns

%set n equal to 1 or 5 or 10
n=1; %horizon (in years) in predictive regressions

TT.ForwardRealizedRP=12*movmean(TT.exret,[0 n*12-1]); %annualized n-year-forward moving realized excess return
ForwardRealizedRP=[TT.ForwardRealizedRP(2:end);NaN];

% Table 9, Panel A or B or C (depending on choice of n above), col. 1
reg=fitlm([a_proxy],ForwardRealizedRP);
tab_data=hac_robust(reg)
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
usedObsBench=usedObs;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

% Table 9, Panel A or B or C (depending on choice of n above), col. 2
reg=fitlm([phi_proxy(usedObsBench).^2],ForwardRealizedRP(usedObsBench));
tab_data=hac_robust(reg)

% Table 9, Panel A or B or C (depending on choice of n above), col. 3
reg=fitlm([a_proxy phi_proxy.^2],ForwardRealizedRP);
tab_data=hac_robust(reg)

% Table 9, Panel A or B or C (depending on choice of n above), col. 4
reg=fitlm([a_proxy phi_proxy phi_proxy.^2],ForwardRealizedRP);
tab_data=hac_robust(reg)