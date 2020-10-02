function plot_pca_loadings(D,varargin)


save_path = '';
for i=1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'SaveDir'
                i = i+1;
                save_path = varargin{i};
        end
    end
end

if numel(D)==2
    D_labels = {'D12';'D3'};
else
    D_labels = {'D1';'D2';'D3'};
end

max_plots = 6;
for i=1:numel(D)
   n_features = cellfun(@(loadings) numel(loadings), D(i).loadings);
   idx = find(n_features>1);
   n_pcs = numel(idx);
   for j=1:n_pcs
       plot_idx = mod(j-1,max_plots)+1;
       if plot_idx==1
           fh=figure;
           fh.Color = [1 1 1];
       end
       ah = subplot(1,max_plots,plot_idx);
       tmp_loadings = D(i).loadings{idx(j)};
       tmp_labels = D(i).loadings_labels{idx(j)};
       [~,p] = sort(tmp_loadings);
       p = fliplr(p');
       tmp_loadings = tmp_loadings(p);
       tmp_labels = tmp_labels(p);
       barh(fliplr(1:numel(p)),tmp_loadings,'EdgeColor','none');
       xlim = get(gca,'XLim');
       xlim = [floor(xlim(1)*10)/10 ceil(xlim(2)*10)/10];
       labels = pretty_labels(fliplr(tmp_labels'));
       set(ah,'XLim',xlim,'YLim',[0 numel(p)+1],'YTick',1:numel(p),...
           'YTickLabels',labels,'fontsize',max_plots,'TickLength',[0 0],...
           'XTick',unique([xlim 0]));
       title(D(i).fields{idx(j)});
       if (plot_idx==max_plots||j==n_pcs) && ~isempty(save_path)
           savefig(fh,sprintf('%s%s_loadings_%i.fig',...
               save_path,D_labels{i},ceil(j/max_plots)));
       end
   end
end