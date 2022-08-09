function [ypred,ysd] = MFGP(f1, f2, x1, x2, y1, y2, Chi, opt)

global ModelInfo

ModelInfo.fe = f2;
ModelInfo.fc = f1;
ModelInfo.Xe = x2;
ModelInfo.Xc = x1;
k = opt.dims;
ModelInfo.ye = y2;
ModelInfo.yc = y1;
ModelInfo.Thetac=fminbnd(@likelihoodc, -3,3);
options = optimoptions(@ga,'display','off');
Params=ga(@likelihoodd, k+1,[],[],[],[],[-6 -10],[6 10],[], options); % [-3 -5],[3 5]
ModelInfo.Thetad=Params(1:k);
ModelInfo.rho=Params (k+1);

buildGP
Xplot = Chi;

ModelInfo.Option='Pred';
pred = ones(length(Xplot(:,end)), 1);
for i=1:length(Xplot(:,end))
    pred(i)=GPpredictor(Xplot(i, :));
end

ModelInfo.Option='RMSE';
standdev = ones(length(Xplot(:,end)), 1);
for i=1:length(Xplot(:,end))
    standdev(i)=GPpredictor(Xplot(i, :));
end

ypred = pred;
ysd = standdev;

end

