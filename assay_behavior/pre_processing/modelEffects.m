function [S,linModels] = modelEffects(S, nuisance_params)

if ~iscell(nuisance_params)
    nuisance_params = {nuisance_params};
end

catvars = {'Plate';'Box';'Tray'};
linModels = cell(numel(cat(1,S.fields)),5);
ct=0;
pct = 0;
pct_update = true;
warning('off','stats:LinearModel:RankDefDesignMat');

for i=1:length(S)
   
    if ~isempty(S(i).data)

        rf = fieldnames(S(i).data);
        rf(strcmp(rf,'filter'))=[];
    
        for j=1:length(rf)

            if mod(pct,24)==0 && pct_update
                figure;
                pct_update = false;
            end
        
            ct=ct+1;
            pf = fieldnames(S(i).meta);
            pf = nuisance_params(any(str_list_contains(nuisance_params,pf),2));
            
            resp_var = S(i).data.(rf{j});
            p = prctile(resp_var,[10 90]);
            filt = resp_var > p(1) & resp_var < p(2);
            
            varnames = [pf; rf(j)];
            vars = cellfun(@(x) S(i).meta.(x)(filt),pf,'UniformOutput',false);
            vars = [vars; resp_var(filt)];
            tbl = table(vars{:},'VariableNames',varnames);
            mdl = fitlm(tbl,'ResponseVar',rf{j},'PredictorVars',...
                pf,'CategoricalVar',catvars(cellfun(@(x) any(strcmp(x,pf)),catvars)));
            sigTerms = mdl.Coefficients{:,4} < 0.01;
            sigTerms(1) = false;

            % remove non-significant terms from model and re-run
            if any(sigTerms)

                sigFields = mdl.CoefficientNames(sigTerms);
                pf = pf(cellfun(@(x) ~isempty(strcmpi(x,sigFields)),pf));
                
                vars = cellfun(@(x) S(i).meta.(x),pf,'UniformOutput',false);
                vars = [vars;S(i).data.(rf{j})];
                varnames = [pf; rf(j)];
                tbl = table(vars{:},'VariableNames',varnames);
                mdl = fitlm(tbl,'ResponseVar',rf{j},'PredictorVars',...
                    pf,'CategoricalVar',catvars(cellfun(@(x) any(strcmp(x,pf)),catvars)));
                
                if ~isempty(pf)

                    label = rf{j};
                    label(label=='_')=' ';
                    pd = makedist('Normal','mu',0,'sigma',0.08);
                    xx = repmat(random(pd,1,numel(vars{1})), 1, 1);
                    opts = {'Marker'; 'o'; 'LineStyle'; 'none';...
                        'MarkerFaceColor'; 'k'; 'MarkerEdgeColor'; 'none';...
                        'MarkerSize'; 3; 'LineWidth'; 1};
                    
                    
                    subplot(4,6,mod(pct,24)+1);
                    pct = pct+1;
                    pct_update = true;
                    if iscell(vars{1})
                        p = cellfun(@(v) find(strcmp(unique(vars{1}),v)), vars{1});
                    else
                        p = vars{1};
                    end
                    plot(p+xx',S(i).data.(rf{j}),opts{:});
                    xlabel(pf);
                    ylabel(label);
                    title(sprintf('%s (%i) - unregressed',S(i).name,S(i).day));
                    set(gca,'XLim',[0 numel(unique(p))+1],...
                        'XTick',unique(p),'XTickLabel',unique(vars{1})); 
                    
                    
                    subplot(4,6,mod(pct,24)+1);
                    pct = pct+1;
                    plot(p+xx',mdl.Residuals.Raw,opts{:});
                    xlabel(pf);
                    ylabel(label);
                    title(sprintf('%s (%i) - regressed',S(i).name,S(i).day));
                    set(gca,'XLim',[0 numel(unique(p))+1],...
                        'XTick',unique(p),'XTickLabel',unique(vars{1})); 
                end

                linModels(ct,1) = {[S(i).name '-' rf{j}]};
                linModels(ct,2) = {mdl.CoefficientNames};
                linModels(ct,3) = {mdl.Coefficients{:,1}};
                linModels(ct,4) = {mdl.Coefficients{:,4}};
                linModels(ct,5) = {S(i).data.(rf{j})};

                S(i).data.(rf{j}) = mdl.Residuals.Raw;
            end

        end
    
    end
    
end
