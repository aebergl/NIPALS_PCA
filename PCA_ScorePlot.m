function fh = PCA_ScorePlot(PCA_Model,XDim,YDim,ColorColumn,CMap,MarkerSymbols,Msize,TitleTxt)
ColorColumnType = 'Group';
LabelPoints=0;
LineWidth = 2;
%ColorColumnType = 'Continues';

if isfield(PCA_Model.Model,'R2X')
    flagPCA = 0;
else
     flagPCA = 1;
end
FontSizeValue = 16;

x = PCA_Model.Model.T(:,XDim);
y = PCA_Model.Model.T(:,YDim);

fh = figure('Name','PCA Score Plot','Color','w','Tag','PCA Score Plot','GraphicsSmoothing','off');
%h=axes(fh,'NextPlot','add','tag','GISTIC bar plot','Box','on','ColorOrder',cmap);

%h = gscatter(x,y,GroupVariable,CMap,MarkerSymbols,Msize);
if strcmp('Group',ColorColumnType)
    GroupVariable = PCA_Model.ObsInfo.Annotation(:,ColorColumn);
%     [~,indx] = sort(GroupVariable);
%     x=x(indx);
%     y=y(indx);
%     GroupVariable=GroupVariable(indx);
    hg = gscatter(x,y,GroupVariable,CMap,MarkerSymbols,Msize);
    ha = hg(1).Parent;
    legend
    [hg.LineWidth] = deal(LineWidth);
    fh.Children(1).AutoUpdate='off';
    fh.Children(1).Title.String=PCA_Model.ObsInfo.AnnotationFields(ColorColumn);
    %fh.Children(1).Title.String = [];

else
    GroupVariable = PCA_Model.ObsInfo.Annotation(:,ColorColumn);
    GroupVariable = str2double(GroupVariable);
    %x = GroupVariable;
    %GroupVariable = ColorColumn;
    hs=scatter(x,y,80,GroupVariable,'Filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor','flat');
    ha = hs.Parent;
    colormap(TurboMap);
    colormap(colorcet('L20'));
    hcb = colorbar;
    hcb.Label.String = PCA_Model.ObsInfo.AnnotationFields(ColorColumn);
    hcb.Label.FontSize = 16;
    ha.Box='on';
    %[CorrValPear,pValPear] =  corr(x,y,'rows','pairwise','type','Pearson')
    %[CorrValSpear,pValSpear] =  corr(x,y,'rows','pairwise','type','Spearman')
end
if flagPCA
    xlabel(sprintf('PC%u (%.1f%%)',XDim,100*PCA_Model.Model.R2(XDim)),'FontSize',FontSizeValue)
    ylabel(sprintf('PC%u (%.1f%%)',YDim,100*PCA_Model.Model.R2(YDim)),'FontSize',FontSizeValue)
else
    %xlabel(sprintf('PC%u (%.1f%%X, %.1f%%Y)',XDim,100*PCA_Model.Model.R2X(XDim),100*PCA_Model.Model.R2Y(XDim)),'FontSize',FontSizeValue)
    %ylabel(sprintf('PC%u (%.1f%%X, %.1f%%Y)',YDim,100*PCA_Model.Model.R2X(YDim),100*PCA_Model.Model.R2Y(YDim)),'FontSize',FontSizeValue)
    xlabel(sprintf('PLS (%u)',XDim),'FontSize',FontSizeValue)
    ylabel(sprintf('PLS (%u)',YDim),'FontSize',FontSizeValue)
end

if LabelPoints
    %MarkerTxt = PCA_Model.ObsInfo.PrimaryId;
    MarkerTxt = PCA_Model.ObsInfo.Annotation(:,2);
    ht = text(x,y,MarkerTxt,'FontSize',10,'VerticalAlignment','Bottom','Interpreter','none');
end
title(TitleTxt)

ha.LineWidth = 1;

min_x = min(x);
max_x = max(x);
min_y = min(y);
max_y = max(y);
nudge_X = range(x) / 30;
nudge_Y = range(y) / 30;

ha.XLim = [min_x - nudge_X, max_x + nudge_X];
ha.YLim= [min_y - nudge_Y, max_y + nudge_Y];

line(ha.XLim,[0 0],'Color','k','LineWidth',1,'LineStyle','-')
line([0 0],ha.YLim,'Color','k','LineWidth',1,'LineStyle','-')
set(gca,'FontSize',16);

if strcmp('Group',ColorColumnType)
    fh.Children(1).Interpreter='none';
%set(findobj(fh.Children(1),'Interpreter','text'),'Interpreter','none')
    legend_objects = findobj(hg,'Type','Line');
    legend_groups={legend_objects(:).DisplayName};
    [~,indx]=sort(legend_groups);
    %indx = [2 4 1 3];
    %indx = 1:numel(legend_groups);
    indx = flip(indx);
   % lh=legend(legend_objects(indx),legend_groups(indx),'location','eastoutside');
    lh=legend(legend_objects,legend_groups,'location','eastoutside');
    %lh=legend(legend_objects(indx),legend_groups(indx));
    %lh=legend(legend_objects,legend_groups);
end
