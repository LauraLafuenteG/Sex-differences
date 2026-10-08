%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Variance-based global sensitivity analysis (Sobol method)  %
%             Implemented by Laura Lafuente-Gracia            %
%                  Last revision: 15/09/2026                  %
% ----------------------------------------------------------- %
%  Main code: Sobol_main.m                                    %
%  Current function: sobol_sampling.m                         %
%    * generates two independent quasi-random sample matrices %
%      that cover the full parameter space                    %
%    * N x D dimension, with N = number or samples            %
%                            D = number of parameters         %
%    * each row is one complete parameter set                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,B] = sobol_sampling(N,LB,UB)

    D = length(LB);     % number of parameters

    % Quasi-Monte Carlo: unlike ordinary Monte Carlo (purely random), Sobol sequences are designed to fill the space  more uniformly              
    % Pure Sobol sequences are deterministic, scrambling introduces a controlled randomization while preserving the space-filling properties
    % So: Sobol sequence + scrambling = less structured, still uniformly distributed, better statistical behaviour

    % Generate one scrambled Sobol sequence in 2D dimensions
    % 1. creates a Sobol sequence generator --> points produced in the unit hypercube [0,1]^2D
    P = sobolset(2*D);
    % 2. scrambling method: randomly permutes digits of the Sobol sequence,
    %                       preserves low discrepancy, avoids some artificial patterns
    P = scramble(P,'MatousekAffineOwen');
    
    % Generate N sample points in the unit hypercube [0,1]^(2D)
    X = net(P,N);
    
    % Split into A and B (size N x D each)
    A = X(:,1:D);
    B = X(:,D+1:2*D);
    
    % Scale from [0,1] to parameter ranges (given by lower and upper bounds)
    LB = LB(:)';
    UB = UB(:)';

    A = LB + A .* (UB-LB);
    B = LB + B .* (UB-LB);

end