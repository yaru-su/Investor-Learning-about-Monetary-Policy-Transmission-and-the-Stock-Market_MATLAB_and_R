function [sumloglik,logLik,stdeps]=LikelihoodFunc_Delta(ConsGrowth,param)

n=length(ConsGrowth);

mudeltabar=abs(param(1,1));
sigmadelta=abs(param(2,1));


for i=1:n-1,   
    eps(i,1)=ConsGrowth(i,1)-(mudeltabar-0.5*sigmadelta^2);
    vareps=(sigmadelta^2);
    stdeps(i,1)=eps(i,1)/sqrt(vareps);
  
    % normal likelihood for whole time period
    epsVec=[eps(i,1)];
    VarMat=[vareps];
    Nv=length(epsVec);

    logLik(i,1)=log(1./((2.*pi)^(Nv/2).*sqrt(det(VarMat))))-1/2*(epsVec*(VarMat\epsVec'));
end

sumloglik=-sum(logLik(2:end,1)); % this the sum of log likelihood (in negative sign since fminsearch minimizes)

end