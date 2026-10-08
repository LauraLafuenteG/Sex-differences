%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Variance-based global sensitivity analysis (Sobol method)  %
%             Implemented by Laura Lafuente-Gracia            %
%                  Last revision: 17/09/2026                  %
% ----------------------------------------------------------- %
%  Sobol convergence analysis:                                %
%    * check convergence to choose N = number of samples      %
%    * creates convergence plots per output                   %
%    * using one nested Sobol sampling                        %
% Additional main code:                                       %
%    * Sobol_main.m: perform global sensitivity analysis      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear;

[p0,LB,UB,parNames,outputNames] = parameters;

D = length(p0);
nOutputs = length(outputNames);

Nlist = [64 128 256 512 1024 2048 4096 8192];
Nmax = max(Nlist);

Results = struct;

tic

% Generate ONE nested Sobol design
[Aall,Ball] = sobol_nestedsampling(Nmax,LB,UB);

% Convergence analysis
for nID = 1:length(Nlist)

    N = Nlist(nID);

    fprintf('\n====================================\n')
    fprintf('N = %d\n',N)
    fprintf('====================================\n')

    % Take nested prefixes
    A = Aall(1:N,:);
    B = Ball(1:N,:);

    YA  = zeros(N,nOutputs);
    YB  = zeros(N,nOutputs);
    YAB = zeros(N,nOutputs,D);

    % Evaluate A
    for i = 1:N
        [~,outVec] = evaluate_model(A(i,:));
        YA(i,:) = outVec;
    end

    % Evaluate B
    for i = 1:N
        [~,outVec] = evaluate_model(B(i,:));
        YB(i,:) = outVec;
    end

    % Evaluate mixed matrices
    for k = 1:D

        fprintf('Parameter %d/%d\n',k,D)

        AB = A;
        AB(:,k) = B(:,k);

        for i = 1:N
            [~,outVec] = evaluate_model(AB(i,:));
            YAB(i,:,k) = outVec;
        end

    end

    % Compute Sobol indices
    Si  = zeros(D,nOutputs);
    STi = zeros(D,nOutputs);

    for outID = 1:nOutputs

        [Si(:,outID),STi(:,outID)] = ...
            sobol_indices( ...
                YA(:,outID), ...
                YB(:,outID), ...
                squeeze(YAB(:,outID,:)));

    end

    Results(nID).N   = N;
    Results(nID).Si  = Si;
    Results(nID).STi = STi;

end

% Save results
save('output/Sobol_convergence.mat','Results','Nlist','parNames','outputNames')


% Quantitative convergence check

Nref = Nlist(end);

fprintf('\n====================================\n')
fprintf('Quantitative convergence check\n')
fprintf('Reference solution: N = %d\n',Nref)
fprintf('====================================\n')

for outID = 1:nOutputs

    STfinal = Results(end).STi(:,outID);

    fprintf('\nOutput: %s\n',outputNames{outID})

    for nID = 1:length(Nlist)-1

        err = abs(Results(nID).STi(:,outID) - STfinal);

        fprintf('N = %5d: max |STi - STi(%d)| = %.5f\n', Nlist(nID), Nref, max(err));

    end
end

% Ranking stability tables
for outID = 1:nOutputs

    filename = sprintf( ...
        'output/Convergence_%s.xlsx', ...
        strrep(outputNames{outID},' ','_'));

    for nID = 1:length(Nlist)

        [~,idx] = sort(Results(nID).STi(:,outID), 'descend');

        T = table( ...
            (1:D)', ...
            parNames(idx), ...
            Results(nID).STi(idx,outID), ...
            'VariableNames',{'Rank','Parameter','STi'});

        writetable(T, filename, 'Sheet',sprintf('N%d',Nlist(nID)));

    end

end

% Convergence plots
for outID = 1:nOutputs

    [~,idxFinal] = sort( Results(end).STi(:,outID), 'descend');

    topIdx = idxFinal(1:min(10,D));

    colors = parula(length(topIdx));

    figure('Position',[100 100 1200 700])
    hold on

    for p = 1:length(topIdx)

        vals = zeros(length(Nlist),1);

        for nID = 1:length(Nlist)
            vals(nID) = Results(nID).STi(topIdx(p),outID);
        end

        plot(Nlist, vals, 'Color',colors(p,:), 'Marker','o', 'LineWidth',1.5);

    end

    xlabel('Base sample size N')
    ylabel('Total-order index ST_i')

    xlim([min(Nlist) max(Nlist)])
    ylim([0 1])

    title(sprintf('ST_i convergence: %s',outputNames{outID}))
    legend(parNames(topIdx),'Location','eastoutside');

    grid on

    exportgraphics(gcf,sprintf('output/Sobol_convergence_%s.png',strrep(outputNames{outID},' ','_')),'Resolution',300);
    close

end

% Finish

elapsedTime = toc;

fprintf('\nFinished convergence analysis.\n');
fprintf('Elapsed time: %.2f seconds\n',elapsedTime);

function [A,B] = sobol_nestedsampling(Nmax,LB,UB)

    D = length(LB);

    % Generate one scrambled Sobol sequence in 2D dimensions.
    %
    % First D columns  -> A
    % Last D columns   -> B
    %
    % Using one sequence allows nested sample sizes:
    % N = 64, 128, ..., Nmax all use prefixes of the same design.

    P = sobolset(2*D);
    P = scramble(P,'MatousekAffineOwen');

    X = net(P,Nmax);

    % Split into A and B
    A = X(:,1:D);
    B = X(:,D+1:2*D);

    % Scale from [0,1] to parameter ranges
    LB = LB(:)';
    UB = UB(:)';

    A = LB + A.*(UB-LB);
    B = LB + B.*(UB-LB);

end