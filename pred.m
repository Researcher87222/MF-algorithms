function f=pred(x)
global ModelInfo

X=ModelInfo.Xc;
y=ModelInfo.yc;
theta=10.^ModelInfo.Thetac;
p=2;

n=size(X,1);

PsicXc=zeros(n,n);
for i=1:n
    for j=i+1:n
        PsicXc(i,j)=exp(-sum(theta.*abs(X(i,:)-X(j,:)).^p));
    end
end
PsicXc=PsicXc+PsicXc'+eye(n)+eye(n).*eps;
[U,p]=chol(PsicXc);
one=ones(n,1);
mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
psi=ones(n,1);

for i=1:n
    psi(i)=exp(-sum(theta.*abs(X(i,:)-x).^p));
end

f=mu+psi'*(U\(U'\(y-one*mu)));