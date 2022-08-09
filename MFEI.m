function [x_new, M, Chi] = MFEI(mu,sigma2,Chi,ObjectiveObservations,opt)

numFidels = opt.numFidels;
FidelityCost = opt.FidelityCost;

alpha1 = zeros(length(Chi(:,end)), numFidels);
alpha2 = zeros(length(Chi(:,end)), numFidels);
alpha3 = zeros(1, numFidels);
sigma2_eps = 0.01; 

best = min(ObjectiveObservations);
EI = compute_ei(best, mu(:,end), sigma2(:,end));

for r = 1:numFidels
    
    co = corr([mu(:,r), mu(:,end)]);
    alpha1(:,r) = co(1,2);
    if r == numFidels
        alpha1(:,r) = 1;
    end
    
    alpha2(:,r) = 1 - sigma2_eps./(sigma2(:,r)+ (sigma2_eps^2));
    alpha3(r) = FidelityCost(end)/FidelityCost(r);
end

MFEI = (EI).*alpha1.*alpha2.*alpha3;



[~, max_index_zero(1,:)] = max(MFEI);
[~, m] = max(max(MFEI));

x_new = Chi(max_index_zero(m), :);
M(1) = m;
Chi(max_index_zero(m),:) = [];

    function ei = compute_ei(best,mu,sigma2)
        sigmas = sqrt(sigma2);
        Z = (best-mu-0.01) ./ sigmas;
        ucdf = normcdf(Z);
        updf = normpdf(Z);
        ei = sigmas .* (Z .* ucdf + updf);
    end



end






