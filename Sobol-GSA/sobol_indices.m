%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Variance-based global sensitivity analysis (Sobol method)  %
%             Implemented by Laura Lafuente-Gracia            %
%                  Last revision: 15/09/2026                  %
% ----------------------------------------------------------- %
%  Main code: Sobol_main.m                                    %
%  Current function: sobol_indices.m                          %
%    * computes Sobol indices for each parameter & output     %
%    * Sobol indices:                                         %
%      * Si = main-effect or first-order index                %
%      * STi = total-effect or total-order index              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Si,STi] = sobol_indices(YA,YB,YAB)
% inputs for one output only, for all parameters k

    D = size(YAB,2);      % number of parameters

    VY = var([YA;YB],1);  % Computes total output variance V(Y) 
                          % using all evaluations from A and B
                          % Creates 2N x 1 samples
                          % Normalizes by N instead of N-1
                 
    Si  = zeros(D,1);     % preallocate
    STi = zeros(D,1);     % preallocate
    
    % compute Sobol indices
    for k = 1:D

        % first-order (main effect)
        Si(k) = mean(YB .* (YAB(:,k)-YA) ) / VY;  % Saltelli's estimator

        % total-order (main effect + interactions)
        STi(k) = mean((YA-YAB(:,k)).^2 ) /(2*VY); % Jansen's estimator
        % Note that only YA and not YB is used because AB = A 
        % that is: AB matrix differs from A only in parameter k (taken from B)
        % So YA-YAB isolates the impact of parameter k: the larger this difference is,
        % the stronger the parameter k influences the output
    end

end
