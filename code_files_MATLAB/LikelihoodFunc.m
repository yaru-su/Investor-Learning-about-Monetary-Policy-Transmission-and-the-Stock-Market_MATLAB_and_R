function [sumloglik,logLik,a,nua,LTMpifinal,stdeps]=LikelihoodFunc(inflation,fedfunds,param)

global Delta

n=length(inflation);

sigmapi=abs(param(1,1));
pibar=param(2,1); 
lambdapi=abs(param(3,1));
sigmaa=abs(param(4,1));
lambdaa=abs(param(5,1));
abar=0;
rNbar=mean(fedfunds);

%initial values
a(1,1)=abar;
nua(1,1)=0.5*sigmaa^2/(2*lambdaa);

for i=1:n-1,   
    %Approximate solution of Ornstein-Uhlenbeck (OU)
    LTMpi(i,1)=pibar-a(i,1)*(fedfunds(i,1)-rNbar);
    eps(i,1)=inflation(i+1,1)-exp(-lambdapi*Delta)*inflation(i,1)-LTMpi(i,1)*(1-exp(-lambdapi*Delta));
    vareps=sigmapi^2/(2*lambdapi)*(1-exp(-2*lambdapi*Delta));
    stdeps(i,1)=eps(i,1)/sqrt(vareps);

    %update a and nua
    vara = ((rNbar-fedfunds(i,1))*lambdapi*nua(i,1)/sigmapi)^2/(2*lambdaa)*(1-exp(-2*lambdaa*Delta)); 
    a(i+1,1)=exp(-lambdaa*Delta)*a(i,1)+abar*(1-exp(-lambdaa*Delta))+sign(rNbar-fedfunds(i,1))*sqrt(vara)*stdeps(i,1);    
       
    nua(i+1,1)=max(0,nua(i,1)+(sigmaa^2-2*lambdaa*nua(i,1)-(rNbar-fedfunds(i,1))^2*lambdapi^2*nua(i,1)^2/sigmapi^2)*Delta);
    
    % normal likelihood for whole time period
    epsVec=[eps(i,1)];
    VarMat=[vareps];
    Nv=length(epsVec);

    logLik(i,1)=log(1./((2.*pi)^(Nv/2).*sqrt(det(VarMat))))-1/2*(epsVec*(VarMat\epsVec'));
end

if mean(LTMpi)+2*std(LTMpi)>quantile(inflation,0.9) || mean(LTMpi)-2*std(LTMpi)<quantile(inflation,0.1) || mean(LTMpi)>mean(inflation)+0.0005 || mean(LTMpi)<mean(inflation)-0.0005 || lambdaa<0.5
    pen=10^10;
else
    pen=0;
end

sumloglik=-sum(logLik(2:end,1))+pen; % this the sum of log likelihood (in negative sign since fminsearch minimizes)

LTMpifinal=pibar-a.*(fedfunds-rNbar);

end