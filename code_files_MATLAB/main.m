clear all; close all; clc

global Delta

Delta=1/12; %data frequency

%% extract data and create timetable

DataGap=readtimetable('data/OutputGapData');%quarterly
DataGap.Properties.VariableNames={'OutputGapQSimple'};
DataGap.Properties.DimensionNames{1}='date';
DataGap.date=DataGap.date+calmonths(2); %downloaded data is for the quarter. Assign the date to Mar, Jun, Sept, Dec instead of Jan, Apr, Jul, Oct
DataGap.OutputGapQ=log(1+DataGap.OutputGapQSimple/100); %annualized continuously compounded output gap
Beg=DataGap.date(1); 
End = datetime(2023, 12, 1);
DataGapMonthly = retime(DataGap(:,{'OutputGapQ'}),Beg:calmonths(1):End,'next');% annualized quarterly output gap to annualized monthly (constant quarterly-value assigned to each month of the quarter) 
DataGapMonthly.Properties.VariableNames={'OutputGap'};

DataFedFunds=readtimetable('data/FEDFUNDS');%monthly
DataFedFunds.Properties.VariableNames={'fedfundsimple'};
DataFedFunds.Properties.DimensionNames{1}='date';
DataFedFunds.fedfund=log(1+DataFedFunds.fedfundsimple/100); %annualized continuously compounded fedfund rate

DataGoyal1=readmatrix("data/GoyalWelch");%monthly
dateStr=num2str(DataGoyal1(:,1));
years=str2num(dateStr(:,1:4));
months=str2num(dateStr(:,5:6));
datetimevector = datetime(years,months,1); %datetime vector format for timetable
DataGoyal2=timetable(datetimevector,DataGoyal1(:,2:end)); %define timetabl
DataGoyal=splitvars(DataGoyal2);
DataGoyal.Properties.DimensionNames{1}='date';%rename columns
DataGoyal.Properties.VariableNames = {'index','d12','e12','bm','tbl','aaa','baa','lty','ntis','rfree','infl','ltr','corpr','svar','csp','sp_vw','sp_vwx'};
DataGoyal.inflation=12*movavg(log(1+DataGoyal.infl),'simple',12,'fill'); %year-over-year continuously compounded inflation (sum 12 consecutive monthly log rates)
DataGoyal.retnom=log(1+DataGoyal.sp_vw); %continuously compounded S&P500 return
DataGoyal=DataGoyal(:,["retnom","d12","index","inflation"]);

DataGDP=readtimetable('data/RealGDP');%annual
DataGDP.Properties.VariableNames={'gdp'};
DataGDP.Properties.DimensionNames{1}='date';
DataGDP.GDPGrowthY=[NaN;log(DataGDP.gdp(2:end,1)./DataGDP.gdp(1:end-1,1))]; %continuously compounded growth rate
DataGDP.date=DataGDP.date+calmonths(11); %assign the date to Dec
Beg=DataGDP.date(1); 
End = datetime(2023, 12, 1);
DataGDPMonthly = retime(DataGDP(:,{'GDPGrowthY'}),Beg:calmonths(1):End,'next');% annual GDP growth to monthly (constant yearly-value assigned to each month of the year)
DataGDPMonthly.Properties.VariableNames={'GDPGrowth'};


%% Download additional variables
DataMartinDaily=readtimetable('data/epbound_Martin.csv');%daily
DataRPMartinMonthly = retime(DataMartinDaily,'monthly','lastvalue'); %last value of the month

DataChabiYoDaily=readtimetable('data/epbound_ChabiYo.csv');%daily
DataRPChabiYoMonthly = retime(DataChabiYoDaily,'monthly','lastvalue'); %last value of the month

DataVIXDaily=readtimetable('data/VIX.csv');%daily
DataVIXMonthly = retime(DataVIXDaily,'monthly','lastvalue'); %last value of the month
DataVIXMonthly.VIX = DataVIXMonthly.VIX/100;

GDPgrowthForecastData = readtable('data/Median_RGDP_Growth.csv');
% Convert to datetime, assuming quarters are Mar, June, Sept, Dec
GDPgrowthForecastData.date = datetime(GDPgrowthForecastData.YEAR, GDPgrowthForecastData.QUARTER * 3, 1);
% Convert to timetable
DataGDPgrowthForecast= table2timetable(GDPgrowthForecastData, 'RowTimes', 'date');
DataGDPgrowthForecast.YEAR=[];
DataGDPgrowthForecast.QUARTER=[];
Beg=min(DataGDPgrowthForecast.date);
End = max(DataGDPgrowthForecast.date);
DataGDPgrowthForecastMonthly = retime(DataGDPgrowthForecast,Beg:calmonths(1):End,'next'); %left-hand side forecasted variable
DataGDPgrowthForecastMonthly=log(1+DataGDPgrowthForecastMonthly./100); %continuously compounded
DataGDPgrowthForecastMonthly=DataGDPgrowthForecastMonthly(:,"DRGDP2");
DataGDPgrowthForecastMonthly = renamevars(DataGDPgrowthForecastMonthly, {'DRGDP2'},{'mudelta'});

%dates at the end have a different format for whatever reason, hence the code below
% Detect structure of the file
opts = detectImportOptions('data/InflationExpectationsClevelandFed.csv');
% Force the first column (date column) to be read as string
opts.VariableTypes{1} = 'string';
% Now read the full table
DataInflationExpectationMonthly = readtable('data/InflationExpectationsClevelandFed.csv', opts);
dates_raw = DataInflationExpectationMonthly.date;
% First attempt: ISO format (e.g., 2022-09-01)
dates1 = datetime(dates_raw, 'InputFormat', 'yyyy-MM-dd', 'Format', 'yyyy-MM-dd', 'Locale', 'en_US');
% Second attempt: short format (e.g., 1-1-2023) for failed cases
bad1 = isnat(dates1);
dates2 = datetime(dates_raw(bad1), 'InputFormat', 'M-d-yyyy', 'Format', 'yyyy-MM-dd', 'Locale', 'en_US');
% Fill in missing values
dates1(bad1) = dates2;
dates = dates1;
DataInflationExpectationMonthly.date=dates;
DataInflationExpectationMonthly=table2timetable(DataInflationExpectationMonthly);

MPUData = readtable('data/MPU.csv'); %monetary policy uncertainty Bloom
MPUData.date = datetime(MPUData.Year, MPUData.Month, 1);
% Convert to timetable
DataMPUMonthly= table2timetable(MPUData, 'RowTimes', 'date');
DataMPUMonthly.Year=[];
DataMPUMonthly.Month=[];
DataMPUMonthly=DataMPUMonthly./100;

CPI5YRForecastData = readtable('data/Median_CPI5YR.csv');
% Convert to datetime, assuming quarters are Mar, June, Sept, Dec
CPI5YRForecastData.date = datetime(CPI5YRForecastData.YEAR, CPI5YRForecastData.QUARTER * 3, 1);
% Convert to timetable
DataCPI5YRForecast= table2timetable(CPI5YRForecastData, 'RowTimes', 'date');
DataCPI5YRForecast.QUARTER=[];
DataCPI5YRForecast.YEAR=[];
Beg=min(DataCPI5YRForecast.date);
End = max(DataCPI5YRForecast.date);
DataCPI5YRForecastMonthly = retime(DataCPI5YRForecast,Beg:calmonths(1):End,'next'); %left-hand side forecasted variable
DataCPI5YRForecastMonthly=DataCPI5YRForecastMonthly./100;


%% merge data in single monthly timetable

TT_list = {
    DataGDPMonthly, DataFedFunds, DataGoyal, DataGapMonthly,DataVIXMonthly, DataRPChabiYoMonthly, ...
    DataRPMartinMonthly, DataGDPgrowthForecastMonthly,DataInflationExpectationMonthly,DataMPUMonthly,DataCPI5YRForecastMonthly
};

% Initialize with the first timetable
TTFull = TT_list{1};
% Loop through the rest and synchronize step by step
for i = 2:length(TT_list)
    TTFull = synchronize(TTFull, TT_list{i});
end

tr=timerange('1954-07-01','2023-12-01','closed');
TT=TTFull(tr,:); %final timetable over '1954-07-01':'2023-12-01' sample
TT.date.Format='yyyy-MM-dd';

%other useful time series
TT.rfnom=TT.fedfund; %nominal risk-free rate
TT.exret=TT.retnom-TT.rfnom/12; %S&P500 excess return
TT.rfreal=TT.rfnom-TT.inflation; %real risk-free rate
TT.dp=log(TT.d12./TT.index); %log dividend price ratio
TT.pd=-TT.dp; %log price-dividend ratio



%% Figure 1

figure(1);
% First panel
subplot(4,1,1);
plot(TT.date, TT.GDPGrowth);
ylabel('GDP growth');
% Second panel
subplot(4,1,2);
plot(TT.date, TT.OutputGap);
ylabel('Output gap');
% Third panel
subplot(4,1,3);
plot(TT.date, TT.inflation);
ylabel('Inflation');
% Fourth panel
subplot(4,1,4);
plot(TT.date, TT.rfnom);
ylabel('Fed funds rate');


%% ML Estimation of (mudeltabar,sigmadelta) mudeltabar not needed in model
options = optimset('TolX',1e-4,'TolFun',1e-4,'Algorithm','interior-point');   
paramdelta0=[0.02;0.02];
[paramdelta,fval]=fminsearch(@(x)LikelihoodFunc_Delta(TT.GDPGrowth,x),paramdelta0,options);
paramdelta;

% Compute t-stat of estimated parameters
[Vdelta]=GetVarMatrixParam_Delta(TT.GDPGrowth,paramdelta);
paramdelta_std = sqrt(diag(Vdelta));
nobs = length(TT.GDPGrowth); 
paramdelta_tstats  = abs(paramdelta./paramdelta_std);
paramdelta_pvalues = 2*(1-tcdf(paramdelta_tstats,nobs-numel(paramdelta)));


%% ML estimation of (sigmay,lambday)
reg=nwest(TT.OutputGap(2:end,1),[ones(length(TT.OutputGap)-1,1) TT.OutputGap(1:end-1,1)],0);
lamby_reg=-log(reg.beta(2))/Delta;

paramy0=[std(TT.OutputGap);lamby_reg];
[paramy,fval]=fminsearch(@(x)LikelihoodFuncy(TT.OutputGap,x),paramy0,options);

% Compute t-stat of estimated parameters
[Vy]=GetVarMatrixParamy(TT.OutputGap,paramy);
paramy_std = sqrt(diag(Vy));
nobs = length(TT.OutputGap); 
paramy_tstats  = abs(paramy./paramy_std);
paramy_pvalues = 2*(1-tcdf(paramy_tstats,nobs-numel(paramy)));


%% ML Estimation of (rNbar,betapi,betay,volr) volr not needed in model
paramr0=[mean(TT.rfnom);0.25;0.15;0.02];
[paramr,fval]=fminsearch(@(x)LikelihoodFunc_Rate(TT.rfnom,TT.inflation,TT.OutputGap,x),paramr0,options);
[sumloglik,logLik,rfTaylor,stdeps]=LikelihoodFunc_Rate(TT.rfnom,TT.inflation,TT.OutputGap,paramr);

% Compute t-stat of estimated parameters
[Vr]=GetVarMatrixParam_Rate(TT.rfnom,TT.inflation,TT.OutputGap,paramr);
paramr_std = sqrt(diag(Vr));
nobs = length(TT.inflation); 
paramr_tstats  = abs(paramr./paramr_std);
paramr_pvalues = 2*(1-tcdf(paramr_tstats,nobs-numel(paramr)));

rNbar_estimated=paramr(1);


%% 5. Inflation and Transmission Mechanism: ML estimation of inflation parameters (sigmapi,pibar,lambdapi,sigmaa,lambdaa)
COV=cov(TT.inflation(2:end,1),TT.inflation(1:end-1,1));
lambdapi0=-1/Delta*log(COV(1,2)/COV(1,1));
sigmapi0=sqrt(var(TT.inflation)*2*lambdapi0);
pibar0=mean(TT.inflation);

param0=[sigmapi0;pibar0;0.1;0.2;0.5];
[param,fval]=fminsearch(@(x)LikelihoodFunc(TT.inflation,TT.fedfund,x),param0,options);
[sumloglik,logLik,a,nua,LTMpi,stdeps]=LikelihoodFunc(TT.inflation,TT.fedfund,param);

% Compute t-stat of estimated parameters
[V]=GetVarMatrixParam(TT.inflation,TT.fedfund,param);
param_std = sqrt(diag(V));
nobs = length(TT.inflation); 
param_tstats  = abs(param./param_std);
param_pvalues = 2*(1-tcdf(param_tstats,nobs-numel(param)));


% Insert ML-Kalman estimated state variables in timetable
TT.a=a;
TT.nua=nua;
TT.LTMpi=LTMpi;
TT.phi=rfTaylor-rNbar_estimated;


%Property of a
change_a=diff(TT.a);
change_pi=diff(TT.inflation);
dummy=(TT.inflation>mean(TT.inflation));
corr(change_a.*dummy(2:end),change_pi.*dummy(2:end));
corr(change_a.*(1-dummy(2:end)),change_pi.*(1-dummy(2:end)));

%Moments of phi
regautocorrphi=nwest(TT.phi(2:end),[ones(length(TT.phi(2:end)),1) TT.phi(1:end-1)],0);
autocorrphi=regautocorrphi.beta(2);
Volphi=std(TT.phi);

%Moments of LTMpi
eg=nwest(LTMpi(2:end,1),[ones(length(LTMpi(2:end,1)),1) LTMpi(1:end-1,1)],0);
autocorLTMpi=reg.beta(2);
volLTMpi=std(LTMpi);

%LTMpi vs. 5-year inflation forecast (footnote 11)
corr(TT.CPI5YR,TT.LTMpi,"rows","complete")
reg=fitlm(TT.LTMpi,TT.CPI5YR);
hac_robust(reg)


%% Figure 2 
figure(2);
% Left y-axis: transmission coefficient
yyaxis left
plot(TT.date, TT.a, 'b'); % blue line
ylabel('Transmission coefficient $\widehat{a}_t$', 'Interpreter', 'latex')

% Right y-axis: inflation
yyaxis right
plot(TT.date, TT.inflation, 'r'); % red line
ylabel('Inflation $\pi_t$', 'Interpreter', 'latex')


%% Figure 3 
figure(3);
% First panel
subplot(2,1,1);
plot(TT.date, TT.phi);
ylabel('Monetary policy stance $\phi_t$', 'Interpreter', 'latex');
% Second panel
subplot(2,1,2);
plot(TT.date, TT.LTMpi);
ylabel('Long-term inflation drift $\widehat{\pi}_t^A$', 'Interpreter', 'latex');



%% Table 3
%sigmadelta,sigmay,lambday,rNbar,betapi,betay,sigmapi,pibar,lambdapi,sigmaa,lambdaa
ParamEstimTable=[paramdelta(2,1) paramdelta_std(2,1);
    paramy paramy_std;
    paramr(1:end-1,1) paramr_std(1:end-1,1);
    param param_std]


%% SAVE ML-KALMAN ESTIMATED STATE VARIABLES
%writetimetable(TT(:,{'GDPGrowth','rfnom','inflation','a','nua','phi','OutputGap','LTMpi'}),'StateVarsOct2024.csv','Delimiter',',')


%% S&P500 real yearly log dividend growth (used to calibrate sigmaD in model)
divgrowth=log(TT.d12(13:end)./TT.d12(1:end-12))-TT.inflation(13:end);
std(divgrowth);

%% Table 4, Data
%real and nominal risk-free rate
meanrfreal=mean(TT.rfreal)
meanrfnom=mean(TT.rfnom)

%S&P500 excess returns
meanSP=12*mean(TT.exret)
volSP=sqrt(12)*std(TT.exret)


%% Table 4, Model is obtained in Mathematica by setting pi=pibar and a=phi=0

%% Model-implied time series
DataModel=readtable('data/ModelImpliedTimeSeriesNov2024'); %these series come from Mathematica
DataModel.date=TT.date;
TTmodel=table2timetable(DataModel);


%% Expected output growth
reg=fitlm([TT.inflation TT.OutputGap],TTmodel.mudelta);
tab_model=hac_robust(reg);
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);
OneStdIncrease_x=std([TT.inflation(usedObs) TT.OutputGap(usedObs)]).*reg.Coefficients.Estimate(2:3)';
OneStd_y=std(TTmodel.mudelta(usedObs));
OneStdIncrease_x/OneStd_y;

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



%% real risk-free rate
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


%% log PD ratio
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




%% market risk premium
%corr between Martin and ChabiYo RPs
corr(TT.RP1, 12*TT.LBR_30, 'Rows', 'complete');

RPmodel=TTmodel.RP;
RPemp=12*TT.LBR_30; %annualized risk premium (1-month-ahead, preference-free) ChabiYo

reg=fitlm([TT.a],RPmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

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


reg=fitlm([TT.phi.^2],RPmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

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



reg=fitlm([TT.a TT.phi.^2],RPmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

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



%% Volatility
Volmodel=TTmodel.vol;
Volemp=TT.VIX;


reg=fitlm([TT.a],Volmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

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


reg=fitlm([TT.phi.^2],Volmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

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


reg=fitlm([TT.a TT.phi.^2],Volmodel);
tab_model=hac_robust(reg);
% Extract the logical index of used observations
usedObs = reg.ObservationInfo.Subset;
% Extract the corresponding dates from TT
datesUsed = TT.Properties.RowTimes(usedObs);
% Get first and last date used in the regression
firstDate = min(datesUsed);
lastDate = max(datesUsed);

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


%% regressions with model-free proxies
phi_proxy=TT.inflation-mean(TT.inflation,'omit');
corr(TT.phi,phi_proxy,'rows','complete')

%compute long-term inflation expectation: mean across horizons 5 to 30
ExpectedInflationMat = [];
for i = 5:1:30
    varname = sprintf('x%dYearExpectedInflation', i);
    ExpectedInflationMat = [ExpectedInflationMat, TT.(varname)];
end
LTExpectedInflationMean=mean(ExpectedInflationMat,2);
LowQuantileInflation=quantile(TT.inflation,0.05);
a_proxy=-LTExpectedInflationMean./max(LowQuantileInflation,TT.inflation);  
corr(TT.a,a_proxy,'rows','complete')

% Table 7, Panel A, col. 1
reg=fitlm(a_proxy,RPemp);
tab_data=hac_robust(reg)
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


%% Impact of nua on RP and Vol
nua_proxy=TT.MPU1; %monetary policy uncertainty, Bloom
IndTime=~isnan(nua_proxy);

%define times with high nua
Times_HighnuaTrue=(TT.nua>quantile(TT.nua(~isnan(TT.nua)),0.5));
Times_Highnua=double(nua_proxy>quantile(nua_proxy(~isnan(nua_proxy)),0.5));
Times_Highnua(isnan(nua_proxy)) = NaN; % Set NaN where nua_proxy is NaN

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


%% Predictive regressions for future realized excess returns

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




