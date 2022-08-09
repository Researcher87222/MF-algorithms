function [NegLnLike,Psi,U] = likelihood(x)


global ModelInfo
X=ModelInfo.Xe;
y=ModelInfo.ye;
theta=10.^x;
n=size(X,1);
one=ones(n,1);
Psi=zeros(n,n);

for i=1:n
    for j=i+1:n
        Psi(i,j)=exp(-sum(theta.*(X(i,:)-X(j,:)).^2));
    end
end

Psi=Psi+Psi'+eye(n)+eye(n).*eps;

[U,p]=chol(Psi);

if p > 0
    NegLnLike=1e4;
else
    
    LnDetPsi=2*sum(log(abs(diag(U))));
    mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
    SigmaSqr=((y-one*mu)'*(U\(U'\ (y-one*mu))))/n;
    NegLnLike=-1*(-(n/2)*log(SigmaSqr) - 0.5*LnDetPsi);
end
end

