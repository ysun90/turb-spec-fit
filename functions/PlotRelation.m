function PlotRelation(d,Para,ParaStr,index)

% PLOTRELATION(D,P,I) plots the relations between
% the specified parameters(P) and all the other parameters
% in the database for the given logical index(I).
% For example:
% PLOTRELATION(d,d.PowerFB,(indOK & indOhmic)).

% get parameter names in the database
ParameterNames = d.Properties.VariableNames';

% plot relations
for ii = 1:length(ParameterNames)
    x = Para(index);
    y = d.(ParameterNames{ii})(index); 
    figure;
    plot(x,y,'b.'); 
    xlabel(ParaStr,'Interpreter','none');
    ylabel(ParameterNames{ii},'Interpreter','none');   
    saveas(gcf,[ParaStr,'_PlotRelation',num2str(ii)],'png');
end