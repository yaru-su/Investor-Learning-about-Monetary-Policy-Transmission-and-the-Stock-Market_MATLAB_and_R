function tab_data=hac_robust(reg)


res = reg.Residuals.Raw(~isnan(reg.Residuals.Raw));
% Number of observations
n = length(res);
% compute autocorrelation and confidence bounds
[acf, lags, bounds] = autocorr(res, 'NumLags', n-1);
% Identify significant lags (where autocorrelation exceeds bounds)
significant = (acf(2:end) > bounds(1)) | (acf(2:end) < -bounds(1));  % skip lag 0
% Find the largest significant lag
max_lag = find(significant, 1, 'last'); 
[EstCoeffCov,se,coeff] = hac(reg,display='off',weights="QS",bandwidth=max_lag);
tstat=coeff./se;
pval=2*(1-tcdf(abs(tstat),reg.NumObservations-reg.NumVariables));

num_rows=length(coeff); %number of rows in table
num_cols=6; %number of columns in table
tab=NaN([num_rows num_cols]); %preallocate with NaNs
tab(:,1)=coeff;
tab(:,2)=se;
tab(:,3)=tstat;
tab(:,4)=pval;
tab(1,5)=reg.Rsquared.Adjusted;
tab(1,6)=reg.NumObservations; 

% Define column headers
var_names = {'Coeff', 'SE', 'tStat', 'pVal', 'AdjR2', 'Nobs'};

% Convert to table
tab_data = array2table(tab, 'VariableNames', var_names);

end