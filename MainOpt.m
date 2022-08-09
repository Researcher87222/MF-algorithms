clear 
close 
clc

fun = struct; 


fun.fid1 = @(x) ((6.*x-2).^2).*sin(12.*x-4); 
fun.fid2 = @(x)  0.5.*(((6.*x-2).^2).*sin(12.*x-4)) + 10.*(x-0.5) -5; 

% Optimization settings
opt = struct; 
opt.dims = 1; % dimension of the design space.
opt.mins = 0; % minimum value for each design variable [min(x1),...,min(xdim)]
opt.maxes = 1; % maximum value for each design variable [max(x1),...,max(xdim)]
opt.max_iters = 2000; % maximum iterations
opt.Budget = 100; % maximum budget for the optimization (100*dim)
opt.grid_size = 10000; % size of the grid to discretize the design space
opt.numFidels = 2; % fidelity levels 
opt.FidelityCost = [0.05, 1]; % fidelity costs  [costLF,..,costHF] 
opt.numSamp = 7; % number of initial design points 
opt.numSampFid = [5,2]; % initial design points for each level of fidelity  [samplesLF,..., samplesHF]
opt.numTests = 1; % number of tests 
opt.AF = 1; % Selection of the multifidelity acquisition function 
% AF = 1 Multifidelity Expected Improvement
% AF = 2 Multifidelity Probability of Improvement

% Run Multifidelity Bayesian Optimization
MultifidBayesOpt(fun,opt)


