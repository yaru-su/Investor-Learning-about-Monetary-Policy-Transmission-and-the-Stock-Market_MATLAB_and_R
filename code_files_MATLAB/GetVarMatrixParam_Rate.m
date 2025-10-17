% The calculation (aproximation) of the first and second derivative of the likelihood
% function is done by a two side finite differences method
%
% Four methods for the calculation of the covariance matrix were implemented here:
% (1) Using second partial derivatives
% (2) Using first partial derivatives (outer product Matrix)
% (3) Using white covariance matrix (NOT USED)
% (4) Using Newey and West Covariance matrix (NOT USED)
%
% References:
%
% HAMILTON, J., D. (1994). Time Series Analysis. Princeton University Press.
%
% NEWEY, B., WEST, D (1987) ‘A Simple, Positive semidefinite, Heteroskedasticity
% and Autocorrelation Consistent Covariance Matrix’ Econometrica, ed. 55, p. 347-370
%
% WHITE, H. (1984) ‘Asymptotic Theory for Econometricians’ New York:
% Academic Press.

function [V]=GetVarMatrixParam_Rate(fedfund,inflation,outputgap,param)
myDelta=1e-15*abs(param);

% First Derivative calculation
[sumloglik,logLikVec1]=LikelihoodFunc_Rate(fedfund,inflation,outputgap,param);
n=length(logLikVec1); 
s=zeros(numel(param),n);

for i=1:numel(param)
  vec=param;
  vec(i)=param(i)+myDelta(i);
 [sumloglik,logLikVec2]=LikelihoodFunc_Rate(fedfund,inflation,outputgap,vec);

 s(i,:)=(logLikVec2-logLikVec1)/(myDelta(i));
end

OP_Matrix=(s*s')/n;

V=inv(OP_Matrix)/n;




end