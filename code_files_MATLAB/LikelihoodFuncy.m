function [sumloglik,logLik,stdeps]=LikelihoodFuncy(outputgap,param)

global Delta

n=length(outputgap);

sigmay=abs(param(1,1));
lambday=param(2,1); 
y=outputgap-mean(outputgap);


for i=1:n-1,   
    
    %Approximate solution of Ornstein-Uhlenbeck (OU)
    eps(i,1)=y(i+1,1)-exp(-lambday*Delta)*y(i,1);
    vareps=sigmay^2/(2*lambday)*(1-exp(-2*lambday*Delta));
    stdeps(i,1)=eps(i,1)/sqrt(vareps);
    
    % normal likelihood for whole time period
    epsVec=[eps(i,1)];
    VarMat=[vareps];
    Nv=length(epsVec);

    logLik(i,1)=log(1./((2.*pi)^(Nv/2).*sqrt(det(VarMat))))-1/2*(epsVec*(VarMat\epsVec'));
end

sumloglik=-sum(logLik(2:end,1)); % this the sum of log likelihood (in negative sign since fminsearch minimizes)

end