clear all; close all; clc

global Delta

Delta=1/12; %data frequency

%% extract data and create timetable

% -----------------------------------
% 1. Output Gap
% ======================================================================
DataGap = readtimetable('data/OutputGapData'); % Read quarterly output gap data
DataGap.Properties.VariableNames = {'OutputGapQSimple'}; % Rename the variable
DataGap.Properties.DimensionNames{1} = 'date'; 
DataGap.date = DataGap.date + calmonths(2); % Shift dates to the middle month of each quarter (Mar, Jun, Sep, Dec)
DataGap.OutputGapQ = log(1 + DataGap.OutputGapQSimple/100); % Convert to annualized continuously compounded output gap

Beg = DataGap.date(1); 
End = datetime(2023, 12, 1);
DataGapMonthly = retime(DataGap(:, {'OutputGapQ'}), Beg:calmonths(1):End, 'next'); % Convert quarterly data to monthly by assigning the quarterly value to each month of the quarter
DataGapMonthly.Properties.VariableNames = {'OutputGap'}; % Rename the variable for monthly timetable

% ======================================================================
% 2. Federal Funds Rate
% ======================================================================
DataFedFunds = readtimetable('data/FEDFUNDS'); % Read monthly federal funds rate data
DataFedFunds.Properties.VariableNames = {'fedfundsimple'}; % Rename the variable
DataFedFunds.Properties.DimensionNames{1} = 'date'; 
DataFedFunds.fedfund = log(1 + DataFedFunds.fedfundsimple/100); % Convert nominal monthly rate to annualized continuously compounded rate

% ======================================================================
% 3. Goyal-Welch monthly data
% ======================================================================
DataGoyal1 = readmatrix("data/GoyalWelch"); % Read raw monthly data from CSV
dateStr = num2str(DataGoyal1(:,1));       % Convert numeric date to string
years = str2num(dateStr(:,1:4));          % Extract year
months = str2num(dateStr(:,5:6));         % Extract month
datetimevector = datetime(years, months, 1); % Create datetime vector for timetable
DataGoyal2 = timetable(datetimevector, DataGoyal1(:,2:end)); % Exclude date column
DataGoyal = splitvars(DataGoyal2);  % Split each variable into separate column
DataGoyal.Properties.DimensionNames{1} = 'date'; % Rename row dimension to 'date'
DataGoyal.Properties.VariableNames = {'index','d12','e12','bm','tbl','aaa','baa','lty','ntis','rfree','infl','ltr','corpr','svar','csp','sp_vw','sp_vwx'};
DataGoyal.inflation = 12*movavg(log(1 + DataGoyal.infl), 'simple', 12, 'fill'); % 12-month log return approximation (annualized)
DataGoyal.retnom = log(1 + DataGoyal.sp_vw); % Compute continuously compounded S&P 500 return
DataGoyal = DataGoyal(:, ["retnom", "d12", "index", "inflation"]); % Keep only relevant columns

% ======================================================================
% 4. Real GDP data
% ======================================================================
DataGDP = readtimetable('data/RealGDP'); % Read annual real GDP data
DataGDP.Properties.VariableNames = {'gdp'}; % Rename the variable to 'gdp'
DataGDP.Properties.DimensionNames{1} = 'date'; % Rename the row dimension to 'date'
DataGDP.GDPGrowthY = [NaN; log(DataGDP.gdp(2:end,1) ./ DataGDP.gdp(1:end-1,1))]; % log(current / previous) gives the continuously compounded growth rate
DataGDP.date = DataGDP.date + calmonths(11); % Assign annual GDP to December
Beg = DataGDP.date(1); % Define the monthly time range
End = datetime(2023, 12, 1);
DataGDPMonthly = retime(DataGDP(:, {'GDPGrowthY'}), Beg:calmonths(1):End, 'next'); % Repeat the annual value for each month in that year
DataGDPMonthly.Properties.VariableNames = {'GDPGrowth'}; % Rename the variable for the monthly timetable


%% Download additional variables

% ======================================================================
% 5. Daily data
% ======================================================================

% Martin Equity Risk Premium (ERP)
DataMartinDaily=readtimetable('data/epbound_Martin.csv');%daily
DataRPMartinMonthly = retime(DataMartinDaily,'monthly','lastvalue'); %last value of the month

% Chabi-Yo Equity Risk Premium (ERP)
DataChabiYoDaily=readtimetable('data/epbound_ChabiYo.csv');%daily
DataRPChabiYoMonthly = retime(DataChabiYoDaily,'monthly','lastvalue'); %last value of the month

% VIX - Market Volatility Index
DataVIXDaily=readtimetable('data/VIX.csv');%daily
DataVIXMonthly = retime(DataVIXDaily,'monthly','lastvalue'); %last value of the month
DataVIXMonthly.VIX = DataVIXMonthly.VIX/100;

% ======================================================================
% 6. GDP growth forecast (median)
% ======================================================================
GDPgrowthForecastData = readtable('data/Median_RGDP_Growth.csv');
% Convert to datetime, assuming quarters are Mar, June, Sept, Dec
GDPgrowthForecastData.date = datetime(GDPgrowthForecastData.YEAR, GDPgrowthForecastData.QUARTER * 3, 1);
% Convert to timetable
DataGDPgrowthForecast= table2timetable(GDPgrowthForecastData, 'RowTimes', 'date');
% Remove YEAR and QUARTER columns since they are no longer needed
DataGDPgrowthForecast.YEAR=[];
DataGDPgrowthForecast.QUARTER=[];
% Define the beginning and end dates for the series
Beg=min(DataGDPgrowthForecast.date);
End = max(DataGDPgrowthForecast.date);
% 'next' assigns each quarterly value to all months of that quarter
DataGDPgrowthForecastMonthly = retime(DataGDPgrowthForecast,Beg:calmonths(1):End,'next'); %left-hand side forecasted variable
% Convert forecast values from percent to continuously compounded rates
DataGDPgrowthForecastMonthly=log(1+DataGDPgrowthForecastMonthly./100); %continuously compounded
% Keep only the relevant column 'DRGDP2' for GDP growth forecast
DataGDPgrowthForecastMonthly=DataGDPgrowthForecastMonthly(:,"DRGDP2");
% Rename the column 'DRGDP2' to 'mudelta' for easier reference
DataGDPgrowthForecastMonthly = renamevars(DataGDPgrowthForecastMonthly, {'DRGDP2'},{'mudelta'});

% ======================================================================
% 7. Inflation Expectations Cleveland Fed
% ======================================================================
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
% Convert to timetable
DataInflationExpectationMonthly=table2timetable(DataInflationExpectationMonthly);

% ======================================================================
% 8. Monetary Policy Uncertainty, MPU
% ======================================================================
MPUData = readtable('data/MPU.csv'); %monetary policy uncertainty Bloom
MPUData.date = datetime(MPUData.Year, MPUData.Month, 1);
% Convert to timetable
DataMPUMonthly= table2timetable(MPUData, 'RowTimes', 'date');
DataMPUMonthly.Year=[];
DataMPUMonthly.Month=[];
% Convert the MPU index from percentage to decimal format
DataMPUMonthly=DataMPUMonthly./100;

% ======================================================================
% 9. Monetary Policy Uncertainty, MPU
% ======================================================================
CPI5YRForecastData = readtable('data/Median_CPI5YR.csv');
% Convert to datetime, assuming quarters are Mar, June, Sept, Dec
CPI5YRForecastData.date = datetime(CPI5YRForecastData.YEAR, CPI5YRForecastData.QUARTER * 3, 1);
% Convert to timetable
DataCPI5YRForecast= table2timetable(CPI5YRForecastData, 'RowTimes', 'date');
% Remove unnecessary columns
DataCPI5YRForecast.QUARTER=[];
DataCPI5YRForecast.YEAR=[];
% Define start and end dates
Beg=min(DataCPI5YRForecast.date);
End = max(DataCPI5YRForecast.date);
% Convert quarterly forecast to monthly by assigning each month's value to the quarterly forecast
DataCPI5YRForecastMonthly = retime(DataCPI5YRForecast,Beg:calmonths(1):End,'next'); %left-hand side forecasted variable
% Convert percent values to decimal
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
