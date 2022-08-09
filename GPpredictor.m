function metric=GPpredictor(x)


global ModelInfo
Xe=ModelInfo.Xe;
Xc=ModelInfo.Xc;
ye=ModelInfo.ye;
yc=ModelInfo.yc;
ne=size(Xe,1);
nc=size(Xc,1);
thetad=10.^ModelInfo.Thetad;
thetac=10.^ModelInfo.Thetac;
p=2;
rho=ModelInfo.rho;
one=ones(nc+ne,1);
y=[yc; ye];
cc=ones(nc,1);
for i=1:nc
    cc(i)=rho*ModelInfo.SigmaSqrc*exp(-sum(thetac.*abs(Xc(i,:)-x).^p));
end
cd=ones(ne,1);
for i=1:ne
    cd(i)=rho^2*ModelInfo.SigmaSqrc*exp(-sum(thetac.*abs(Xe(i,:)-x).^p))+ModelInfo.SigmaSqrd*exp(-sum(thetad.*abs(Xe(i,:)-x).^p));
end
c=[cc;cd];
f=ModelInfo.mu+c'*(ModelInfo.UC\(ModelInfo.UC'\(y-one.*ModelInfo.mu)));
if strcmp(ModelInfo.Option,'Pred')==0
    SSqr=rho^2*ModelInfo.SigmaSqrc+ModelInfo.SigmaSqrd-c'*(ModelInfo.UC\(ModelInfo.UC'\c));
    s=abs(SSqr)^0.5;
    if strcmp(ModelInfo.Option,'RMSE')==0
        yBest=min(ye);
        if strcmp(ModelInfo.Option,'NegProbImp')==1
            ProbImp=0.5+0.5*erf((1/sqrt(2))*((yBest-f)/s));
        else
            EITermOne=(yBest-f)*(0.5+0.5*erf((1/sqrt(2))*((yBest-f)/s)));
            EITermTwo=s*(1/sqrt(2*pi))*exp(-(1/2)*((yBest-f)^2/SSqr));
            ExpImp=log10(EITermOne+EITermTwo+realmin);
        end
    end
end
if strcmp(ModelInfo.Option,'Pred')==1
    metric=f;
elseif strcmp(ModelInfo.Option,'RMSE')==1
    metric=s;
elseif strcmp(ModelInfo.Option,'NegLogExpImp')==1
    metric=-ExpImp;
elseif strcmp(ModelInfo.Option,'NegProbImp')==1
    metric=-ProbImp;
end
