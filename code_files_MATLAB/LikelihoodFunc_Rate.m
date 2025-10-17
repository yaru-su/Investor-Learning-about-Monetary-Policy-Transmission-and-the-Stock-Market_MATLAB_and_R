function [sumloglik,logLik,rfTaylor,stdeps]=LikelihoodFunc_Rate(fedfund,inflation,outputgap,param)

n=length(fedfund);

rNbar=param(1,1);
betapi=abs(param(2,1));
betay=abs(param(3,1)); 
volr=abs(param(4,1)); %vol of residual. Not needed in the model

meaninflation=mean(inflation);
meanoutputgap=mean(outputgap);

for i=1:n,   
    rfTaylor(i,1)=rTaylor(rNbar,betapi,betay,inflation(i,1),outputgap(i,1),meaninflation,meanoutputgap);
    eps(i,1)=fedfund(i,1)-rfTaylor(i,1);
    vareps=volr^2;
    stdeps(i,1)=eps(i,1)/sqrt(vareps);
  
    % normal likelihood for whole time period
    epsVec=[eps(i,1)];
    VarMat=[vareps];
    Nv=length(epsVec);

    logLik(i,1)=log(1./((2.*pi)^(Nv/2).*sqrt(det(VarMat))))-1/2*(epsVec*(VarMat\epsVec'));
end


if abs(mean(rfTaylor)-mean(fedfund))>0.00005 || std(betapi*inflation)>4*std(betay*outputgap)
    pen=10^10;
else
    pen=0;
end
sumloglik=-sum(logLik(2:end,1))+pen; % this the sum of log likelihood (in negative sign since fminsearch minimizes)


end