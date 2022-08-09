function NegLnLikec=likelihoodc(x)



global ModelInfo
Xc=ModelInfo.Xc;
yc=ModelInfo.yc;
nc=size(Xc,1);
thetac=10.^x;
p=2;
one=ones(nc,1);
PsicXc=zeros(nc,nc);
for i=1:nc
    for j=i+1:nc
        PsicXc(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xc(j,:)).^p));
    end
end
PsicXc=PsicXc+PsicXc'+eye(nc)+eye(nc).*eps;
[U,p]=chol(PsicXc);
if p>0
    NegLnLikec=100;
else
    LnDetPsicXc=2*sum(log(abs(diag(U))));
    muc=(one'*(U\(U'\yc)))/(one'*(U\(U'\one)));
    SigmaSqrc=(yc-one.*muc)'*(U\(U'\(yc-one.*muc)))/nc;
    NegLnLikec=-1*(-(nc/2)*log(SigmaSqrc)-0.5*LnDetPsicXc);
end

