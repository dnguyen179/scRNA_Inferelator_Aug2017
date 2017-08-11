function histScatterLog(xstuff,ystuff,nbins,stepSize,label,fontSize,figBase)
% use hist3d to get a 2D scatter density plot
% Inputs:
%   xstuff -- column vector of x values
%   ystuff -- column vector of y values
%   nbins -- number of bins for the 2D histogram (will be symmetric in x and
%       y)
%   stepSize -- how often to add an xTick mark or yTick mark on the plot
%   label -- Matlab object
%       label.x = string for the x label
%       label.y = string for the y label
%       label.tit = string for the title
%   fontSize -- fontSize for the figure
%   figBase -- (optional) a filename for saving the figure as .fig + .pdf
% Outputs:
%   a plot and, optionally, saved versions of that plot


%% debugging:
% xstuff=meanIntens;
% ystuff=dataCv;
% nbins = 101;
% stepSize = 10;
% label.x = '\mu_{intens}';
% label.y = 'CV';
% label.tit = 'CV vs. \mu_{intens}';
% fontSize = 14;

%%

totDupeProbs = length(xstuff);
tickInds = [1:stepSize:nbins];
% because of the flipud operation, have to make sure that nbins is
% % symmetric
% if max(tickInds) < nbins
%     tickInds = [tickInds nbins];
% end
% nbins = nbins + 1;
% figure(3), clf
[values, N] = hist3([xstuff ystuff],[nbins nbins]);
imagesc(flipud(log10(values'/totDupeProbs)))
%     imagesc(flipud(values'/totDupeProbs))
%     colormap(flipud(gray))
colormap(flipud(gray))
set(gca,'FontSize',fontSize)
%         colorbar('YTick',
c = colorbar;
c.Label.String = 'log_{10}(Frequency)';
xlabel(label.x,'FontSize',fontSize)
ylabel(label.y,'FontSize',fontSize)

% set(gca,'XTick',tickInds,'XTickLabel',num2str(N{1}(tickInds)','%1.1f'),'FontSize',fontSize)
set(gca,'XTick',tickInds,'XTickLabel',roundstring1((N{1}(tickInds)')),'FontSize',fontSize)
% xticklabel_rotate(tickInds,90,cellstr(roundstring1((N{1}(tickInds)'))),'FontSize',fontSize)
set(gca,'YTick',sort(nbins-tickInds),'YTickLabel',num2str(flipud(N{2}(tickInds)'),'%1.2f'),'FontSize',fontSize)
%         set(gca,'YTick',tickInds,'YTickLabel',num2str(10.^(flipud(N{2}(tickInds)')),'%1.0f'),'FontSize',fontSize)
% hold on
% ax = axis();
% [minVal minInd] = min((10.^N{1} - maxProbIntens).^2);
% plot(minInd*[1 1],[ax(3) ax(4)],'k:','LineWidth',3)
% [minVal minInd] = min((flipud(N{2}') - medCorrCut).^2);
% plot([ax(1) ax(2)],minInd*[1 1],'k:','LineWidth',3)

title([label.tit],...
    'FontSize',fontSize+2)

if figBase
    figName = [figBase];
    save2pdf([figName '.pdf'],gcf,150)
    saveas(gcf,figName,'fig')
    disp([figName '.pdf + .fig'])
end