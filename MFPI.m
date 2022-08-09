function [x_new, M, Chi] = MFPI(mu,sigma2,Chi,ObjectiveObservations,X1,X2,opt)

global ModelInfo
numFidels = opt.numFidels;
FidelityCost = opt.FidelityCost;
thetad=10.^ModelInfo.Thetad;
thetac=10.^ModelInfo.Thetac;
Rho=ModelInfo.rho;
eta=zeros(length(Chi(:,end)),numFidels);
fmin = min(ObjectiveObservations);
sigmas=sqrt(sigma2);
PI = normcdf((fmin-mu(:,end))./sigmas(:,end));
R=ones(length(mu(:,end)),numFidels);

for j= 1: numFidels
    
    C= corr([mu(:,j), mu(:,end)]);
    cor(j)=C(1,2);
    
    if j==numFidels
        cor(j)=1;
    end
    
    CR(j) = FidelityCost(end)/FidelityCost(j);
    
end


for j=1:length(Chi(:,end))
    for k=1: length(X1(:,end))
        R(j,k)=1.-exp(-sum(thetac.*(Chi(j,:)-X1(k,:)).^2));
    end
    eta(j,1) = prod(R(j,:));
end

thetae= Rho*thetac+thetad;

for j=1:length(Chi(:,end))
    for k=1: length(X2(:,end))
        R(j,k)=1.-exp(-sum(thetae.*(Chi(j,:)-X2(k,:)).^2));
    end
    eta(j,2) = prod(R(j,:));
end

MFPI = PI.*cor.*CR.*eta;




[~, max_index_zero(1,:)] = max(MFPI);
[~, m] = max(max(MFPI));

x_new = Chi(max_index_zero(m), :);
M(1) = m;
Chi(max_index_zero(m),:) = [];

end


