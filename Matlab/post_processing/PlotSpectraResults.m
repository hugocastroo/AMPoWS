%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the an estimation of the REWS - Estimate
%from the KF in the ROSCO Controller
%------------------------------------------------------------
%V1.0 2022.11.12 - HC
% ----------------------------------    
function PlotSpectraResults(FrequenciesToPlot,SignalsToPlot,SignalNames,SavePlot,URef,xText,yText,tText,logViewx,logViewy)

    %plot the results
    figure
    hold on;grid on;box on
    xlabel(xText,'FontSize', 20)
    ylabel(yText,'FontSize', 20)
    title(tText);

    if logViewx
        set(gca,'xScale','log')
    end

    if logViewy
        set(gca,'yScale','log')
    end

    for i = 1:length(FrequenciesToPlot)
        plot(FrequenciesToPlot{i},SignalsToPlot{i},'Linewidth',2)
    end
    
    LegendCell = SignalNames;
    legend(LegendCell)

    %Save plot as pdf for further use if desired
    if SavePlot
        pdf = gcf;
        set(pdf,'PaperOrientation','landscape');
        set(pdf,'PaperUnits','normalized');
        set(pdf,'PaperPosition', [0 0 1 1]);
        pdfpath = ['Windspeed_', num2str(URef), 'ms'];
        saveas(gcf,pdfpath,'pdf');
    end