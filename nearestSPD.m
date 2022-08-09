function Ahat = nearestSPD(A)

[r,c] = size(A);
if r ~= c
  error('A not square matrix.')
elseif (r == 1) && (A <= 0)
  Ahat = eps;
  return
end
B = (A + A')/2;
[U,Sigma,V] = svd(B);
H = V*Sigma*V';
Ahat = (B+H)/2;
Ahat = (Ahat + Ahat')/2;
p = 1;
k = 0;
while p ~= 0
  [R,p] = chol(Ahat);
  k = k + 1;
  if p ~= 0
    mineig = min(eig(Ahat));
    Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
  end
end






