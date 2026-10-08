%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Variance-based global sensitivity analysis (Sobol method)  %
%             Implemented by Laura Lafuente-Gracia            %
%                  Last revision: 15/09/2026                  %
% ----------------------------------------------------------- %
%  Main code: Sobol_main.m                                    %
%  Current function: plot_sobol.m                             %
%    * plot main- and total-effect indexes per output         %
%    * plot interactions per output                           %
%    * create heatmap of total-effect index per output        %
%    * create heatmap of interactions per output              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_sobol(Si,STi,parNames,outputNames,N)

    nOutputs = size(Si,2);
    
    % % plot main-effect and total-effect results per output
    % for outID = 1:nOutputs
    % 
    %     [~,idx] = sort(STi(:,outID),'descend');
    % 
    %     figure('Position',[100 100 1400 600])
    % 
    %     bar([Si(idx,outID) STi(idx,outID)])
    % 
    %     xticks(1:length(parNames))
    %     xticklabels(parNames(idx))
    %     xtickangle(90)
    %     ylim([0 1])
    % 
    %     ylabel('Sobol index')
    %     legend('First-order','Total-order')
    %     title(sprintf('Sobol indices: %s (N = %d)', outputNames{outID}, N))
    % 
    %     exportgraphics(gcf, sprintf('output/Sobol_%s_N%d.png', strrep(outputNames{outID},' ','_'),N), 'Resolution',300)
    %     close
    % 
    % end

    % plot with 3 bars
    for outID = 1:nOutputs
            
        SiPlot = max(Si(:,outID),0);
        InteractionPlot = STi(:,outID) - SiPlot;
        
        [~,idx] = sort(STi(:,outID),'descend');
        figure('Position',[100 100 1400 600])
        
        bar([SiPlot(idx), STi(idx,outID), InteractionPlot(idx)])

        xticks(1:length(parNames))
        xticklabels(parNames(idx))
        xtickangle(90)
        ylim([0 1])
    
        ylabel('Sobol index')
        legend('First-order','Total-order','Interactions')
        title(sprintf('Sobol indices: %s (N = %d)', outputNames{outID}, N))
    
        exportgraphics(gcf, sprintf('output/Sobol3bars_%s_N%d.png', strrep(outputNames{outID},' ','_'),N), 'Resolution',300)
        close
    
    end


    % % plot with stacked bars
    % for outID = 1:nOutputs
    % 
    %     Interaction = STi(:,outID)-Si(:,outID);
    %     [~,idx] = sort(STi(:,outID),'descend');
    % 
    %     figure('Position',[100 100 1400 600])
    % 
    %     bar([Si(idx,outID) Interaction(idx)], 'stacked')
    % 
    %     xticks(1:length(parNames))
    %     xticklabels(parNames(idx))
    %     xtickangle(90)
    %     ylim([0 1])
    % 
    %     ylabel('Total-order Sobol index')
    %     legend('Main effect','Interactions')
    % 
    %     title(sprintf('Sobol decomposition: %s', outputNames{outID}))
    % 
    %     exportgraphics(gcf, sprintf('output/Sobol_stacked_%s_N%d.png', strrep(outputNames{outID},' ','_'),N), 'Resolution',300)
    %     close
    % 
    % end


    % % plot interaction results per output
    % for outID = 1:nOutputs
    % 
    %     Interaction = STi(:,outID)-Si(:,outID);
    %     [~,idx] = sort(Interaction,'descend');
    % 
    %     figure
    % 
    %     bar(Interaction(idx))
    % 
    %     xticks(1:length(parNames))
    %     xticklabels(parNames(idx))
    %     xtickangle(90)
    % 
    %     ylabel('ST_i - S_i')
    % 
    %     title(sprintf('Parameter interactions: %s (N = %d)', outputNames{outID}, N))
    % 
    %     exportgraphics(gcf, sprintf('output/interactions_%s_N%d.png', strrep(outputNames{outID},' ','_'),N), 'Resolution',300)
    % 
    %     close
    % 
    % end

    %% Heatmaps
    peakIdx = 1:8;
    ttpIdx  = 9:15;
    InteractionFraction = (STi-Si)./max(STi,1e-12);

    % Order parameters & save Sobol ordering for Matlab and for R
    MeanSTi = mean(STi,2);
    [~,idx] = sort(MeanSTi,'descend');

    sortedParNames = parNames(idx); 
    save('output/Sobol_sorting.mat','idx','sortedParNames')
    RankingTable = table((1:length(idx))',sortedParNames(:),MeanSTi(idx),'VariableNames',{'Rank','Parameter','MeanSTi'});
    writetable(RankingTable,sprintf('output/Sobol_ranking_N%d.xlsx',N));

    % Total-effect heatmaps: all outputs
    figure('WindowState','maximized')
    imagesc(STi(idx,:)')
    colormap(parula)
    colorbar
    clim([0 1])
    xticks(1:length(parNames))
    xticklabels(parNames(idx))
    xtickangle(90)
    yticks(1:length(outputNames))
    yticklabels(outputNames)
    xlabel('Parameter')
    ylabel('Output')
    title(sprintf('Sobol total-order indices (ST_i) (N = %d)',N))
    
    exportgraphics(gcf,sprintf('output/heatmap_total_N%d.png',N),'Resolution',300)
    close

    % export data for R
    T = array2table(STi(idx,:)','VariableNames', matlab.lang.makeValidName(parNames(idx)));
    T = addvars(T,outputNames(:),'Before',1,'NewVariableNames',"Output");
    writetable(T, sprintf('output/heatmap_total_all_N%d.xlsx',N));

    % Total-effect heatmaps: peak outputs
    figure('WindowState','maximized')
    imagesc(STi(idx,peakIdx)')
    colormap(parula)
    colorbar
    clim([0 1])
    xticks(1:length(parNames))
    xticklabels(parNames(idx))
    xtickangle(90)
    yticks(1:length(peakIdx))
    yticklabels(outputNames(peakIdx))
    xlabel('Parameter')
    ylabel('Peak output')
    title(sprintf('Sobol total-order indices (ST_i): peaks (N = %d)',N))

    exportgraphics(gcf,sprintf('output/heatmap_total_peaks_N%d.png',N),'Resolution',300)
    close

    % export data for R
    T = array2table(STi(idx,peakIdx)','VariableNames', matlab.lang.makeValidName(parNames(idx)));
    T = addvars(T,outputNames(peakIdx),'Before',1,'NewVariableNames',"Output");
    writetable(T,sprintf('output/heatmap_total_peaks_N%d.xlsx',N));

    % Total-effect heatmaps: time-to-peak outputs
    figure('WindowState','maximized')
    imagesc(STi(idx,ttpIdx)')
    colormap(parula)
    colorbar
    clim([0 1])
    xticks(1:length(parNames))
    xticklabels(parNames(idx))
    xtickangle(90)
    yticks(1:length(ttpIdx))
    yticklabels(outputNames(ttpIdx))
    xlabel('Parameter')
    ylabel('Time-to-peak output')
    title(sprintf('Sobol total-order indices (ST_i): time-to-peak (N = %d)',N))

    exportgraphics(gcf,sprintf('output/heatmap_total_ttp_N%d.png',N),'Resolution',300)
    close

    % export data for R
    T = array2table(STi(idx,ttpIdx)','VariableNames', matlab.lang.makeValidName(parNames(idx)));
    T = addvars(T,outputNames(ttpIdx),'Before',1,'NewVariableNames',"Output");
    writetable(T,sprintf('output/heatmap_total_ttp_N%d.xlsx',N));

    % Interaction-fraction heatmaps: all outputs    
    figure('WindowState','maximized')
    imagesc(InteractionFraction(idx,:)')
    colormap(parula)
    colorbar
    clim([0 1])
    xticks(1:length(parNames))
    xticklabels(parNames(idx))
    xtickangle(90)
    yticks(1:length(outputNames))
    yticklabels(outputNames)
    xlabel('Parameter')
    title(sprintf('Interaction fraction (N = %d)',N))
    
    exportgraphics(gcf,sprintf('output/heatmap_interactions_N%d.png',N),'Resolution',300)
    close

    % export data for R
    T = array2table(InteractionFraction(idx,:)','VariableNames', matlab.lang.makeValidName(parNames(idx)));
    T = addvars(T,outputNames(:),'Before',1,'NewVariableNames',"Output");
    writetable(T,sprintf('output/heatmap_interactions_all_N%d.xlsx',N));

    % Interaction-fraction heatmaps: peak outputs
    figure('WindowState','maximized')
    imagesc(InteractionFraction(idx,peakIdx)')
    colormap(parula)
    colorbar
    clim([0 1])
    xticks(1:length(parNames))
    xticklabels(parNames(idx))
    xtickangle(90)
    yticks(1:length(peakIdx))
    yticklabels(outputNames(peakIdx))
    xlabel('Parameter')
    ylabel('Peak output')
    title(sprintf('Interaction fraction: peaks (N = %d)',N))

    exportgraphics(gcf,sprintf('output/heatmap_interactions_peaks_N%d.png',N),'Resolution',300)
    close

    % export data for R
    T = array2table(InteractionFraction(idx,peakIdx)','VariableNames', matlab.lang.makeValidName(parNames(idx)));
    T = addvars(T,outputNames(peakIdx),'Before',1,'NewVariableNames',"Output");
    writetable(T,sprintf('output/heatmap_interactions_peaks_N%d.xlsx',N));

    % Interaction-fraction heatmaps: time-to-peak outputs
    figure('WindowState','maximized')
    imagesc(InteractionFraction(idx,ttpIdx)')
    colormap(parula)
    colorbar
    clim([0 1])
    xticks(1:length(parNames))
    xticklabels(parNames(idx))
    xtickangle(90)
    yticks(1:length(ttpIdx))
    yticklabels(outputNames(ttpIdx))
    xlabel('Parameter')
    ylabel('Time-to-peak output')
    title(sprintf('Interaction fraction: time-to-peak (N = %d)',N))

    exportgraphics(gcf,sprintf('output/heatmap_interactions_ttp_N%d.png',N),'Resolution',300)
    close

    % export data for R
    T = array2table(InteractionFraction(idx,ttpIdx)','VariableNames', matlab.lang.makeValidName(parNames(idx)));
    T = addvars(T,outputNames(ttpIdx),'Before',1,'NewVariableNames',"Output");
    writetable(T,sprintf('output/heatmap_interactions_ttp_N%d.xlsx',N));

    %% Interaction-fraction heatmap: sex-specific parameters only
    sexParNames = { ...
        'ds',...
        'k01',...
        'k02',...
        'k1',...
        'k2',...
        'k3',...
        'ke1',...
        'ke2',...
        'klb',...
        'kls',...
        'kpb',...
        'kps',...
        'Mmax'};
    
    sexIdx = find(ismember(parNames,sexParNames));
    sexIdxOrdered = idx(ismember(idx,sexIdx));

    figure('WindowState','maximized')
    imagesc(InteractionFraction(sexIdxOrdered,:)')
    
    colormap(parula)
    colorbar
    clim([0 1])
    
    xticks(1:length(sexIdxOrdered))
    xticklabels(parNames(sexIdxOrdered))
    xtickangle(90)
    
    yticks(1:length(outputNames))
    yticklabels(outputNames)
    
    xlabel('Sex-specific parameter')
    ylabel('Output')
    title(sprintf('Interaction contribution of the 13 sex-specific parameters (N = %d)',N))
    
    exportgraphics(gcf,sprintf('output/heatmap_interactions_sexspecific_N%d.png',N),'Resolution',300)
    close

end