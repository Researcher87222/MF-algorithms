function NegLnLiked = likelihoodd(x)


global ModelInfo
Xe=ModelInfo.Xe;
Xc=ModelInfo.Xc;
ye=ModelInfo.ye;
yc=ModelInfo.yc;
[ne,k]=size(Xe);
thetad=10.^x(1:k);
rho=x(k+1);
one=ones(ne,1);

PsidXe=zeros(ne,ne);

for i=1:ne
    for j=i+1:ne
        PsidXe(i,j)=exp(-sum(thetad.*(Xe(i,:) -Xe(j,:)).^2));
    end
end

PsidXe=PsidXe+PsidXe'+eye(ne)+eye(ne).*eps;

[U,p]=chol(PsidXe);

if p > 0
    NegLnLiked=1e4;
else
    
    LnDetPsidXe=2*sum(log(abs(diag(U))));
    
    if length(yc) - ne +1 > 0
        d = ye-rho.*yc(end-ne+1:end);
    else
        for j = 1:length(Xe(:,1))
            yc(j,1) = pred(Xe(j,:));
        end
        d = ye-rho.*yc;
    end
    
    mud=(one'*(U\(U'\d)))/(one'*(U\(U'\one)));
    SigmaSqrd=(d-one.*mud)'*(U\(U'\(d-one.*mud)))/ne;
    NegLnLiked=-1*(-(ne/2)*log(SigmaSqrd)-0.5*LnDetPsidXe);
end
end

