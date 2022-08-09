function buildGP

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

one=ones(ne+nc,1);
y=[yc; ye];
PsicXc=zeros(nc,nc);
for i=1:nc
    for j=i+1:nc
        PsicXc(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xc(j,:)).^p));
    end
end
ModelInfo.PsicXc=PsicXc+PsicXc'+eye(nc)+eye(nc).*eps;
ModelInfo.UPsicXc=chol(nearestSPD(ModelInfo.PsicXc)); % changed

PsicXe=zeros(ne,ne);
for i=1:ne
    for j=i+1:ne
        PsicXe(i,j)=exp(-sum(thetac.*abs(Xe(i,:)-Xe(j,:)).^p));
    end
end
ModelInfo.PsicXe=PsicXe+PsicXe'+eye(ne)+eye(ne).*eps;
ModelInfo.UPsicXe=chol(nearestSPD(ModelInfo.PsicXe)); % changed


PsicXcXe=zeros(nc,ne);
for i=1:nc
    for j=1:ne
        PsicXcXe(i,j)=exp(-sum(thetac.*abs(Xc(i,:)-Xe(j,:)).^p));
    end
end
if nc-ne>0
    ModelInfo.PsicXcXe=PsicXcXe+[zeros(nc-ne,ne);eye(ne)].*eps;
else
    ModelInfo.PsicXcXe=PsicXcXe + eps;
end

ModelInfo.PsicXeXc=ModelInfo.PsicXcXe';


PsidXe=zeros(ne,ne);
for i=1:ne
    for j=i+1:ne
        PsidXe(i,j)=exp(-sum(thetad.*abs(Xe(i,:)-Xe(j,:)).^p));
    end
end
ModelInfo.PsidXe=PsidXe+PsidXe'+eye(ne)+eye(ne).*eps;
ModelInfo.UPsidXe=chol(nearestSPD(ModelInfo.PsidXe)); %changed


ModelInfo.muc=(ones(nc,1)'*(ModelInfo.UPsicXc\(ModelInfo.UPsicXc'\yc)))/(ones(nc,1)'*(ModelInfo.UPsicXc\(ModelInfo.UPsicXc'\ones(nc,1))));

if length(yc) - ne +1 > 0
    ModelInfo.d = ye-rho.*yc(end-ne+1:end);
else
    for j = 1:length(Xe(:,1))
        ycc(j,1) = pred(Xe(j,:));
    end
    ModelInfo.d = ye-rho.*ycc;
end

ModelInfo.mud=(ones(ne,1)'*(ModelInfo.UPsidXe\(ModelInfo.UPsidXe'\ModelInfo.d)))/(ones(ne,1)'*(ModelInfo.UPsidXe\(ModelInfo.UPsidXe'\ones(ne,1))));

ModelInfo.SigmaSqrc=(yc-ones(nc,1).*ModelInfo.muc)'*(ModelInfo.UPsicXc\(ModelInfo.UPsicXc'\(yc-ones(nc,1).*ModelInfo.muc)))/nc;
ModelInfo.SigmaSqrd=(ModelInfo.d-ones(ne,1).*ModelInfo.mud)'*(ModelInfo.UPsidXe\(ModelInfo.UPsidXe'\(ModelInfo.d-ones(ne,1).*ModelInfo.mud)))/ne;

ModelInfo.C=[ModelInfo.SigmaSqrc*ModelInfo.PsicXc rho*ModelInfo.SigmaSqrc*ModelInfo.PsicXcXe;
    rho*ModelInfo.SigmaSqrc*ModelInfo.PsicXeXc rho^2*ModelInfo.SigmaSqrc*ModelInfo.PsicXe+ModelInfo.SigmaSqrd*ModelInfo.PsidXe];
ModelInfo.UC=chol(nearestSPD(ModelInfo.C));

ModelInfo.mu=(one'*(ModelInfo.UC\(ModelInfo.UC'\y)))/(one'*(ModelInfo.UC\(ModelInfo.UC'\one)));

