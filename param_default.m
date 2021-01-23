function Param = param_default

% Default prameters for the example femoral FBI dataset. May need to adjust
% according to data sparsity etc.

% Iteration number
Param.nInner=4;
Param.nBreg=3;

% Used for k-space-subtraction methods
Param.sub.lambda=1e-4;  % Increase if data sparsity is low.
Param.sub.threshTV=2e-1;
Param.sub.nu=1e-5;

% Used for magnitude-subtraction methods
Param.mag.lambda=5e-3;  
Param.mag.threshTV=6e-2;
Param.mag.nu=1e-5;

Param.IC.alpha=1;   % >0. Sensitivity to background tissues with large signal intensity. Reduce if background suppression is too much.

end

