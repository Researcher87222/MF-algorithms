function [L] = stableCholesky(K)

  try
    L = chol(K, 'lower');
    return;
  catch err
  end
  

  diagPower = min( ceil(log10(abs(min(diag(K)))))-1, -11);
  if ~(abs(diagPower) < inf)
    diagPower = -10;
  end

  success = false;
  K = K + (10^diagPower) * eye(size(K));
  while ~success
    try
      L = chol(K, 'lower');
      success = true;
    catch err
      diagPower = diagPower + 1; 
      K = K + (10^diagPower) * eye(size(K));
    end
  end

end
