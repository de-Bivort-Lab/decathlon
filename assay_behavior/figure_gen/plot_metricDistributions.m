function plot_metricDistributions(D, varargin)
% plot metric kde for decathlon data struct or array of decathlon data structs D
% If D is an array of decathlon data structs, metrics will be paired
% across data structs and plotted on the same axis

% parse inputs
labels = '';
for i=1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'Labels'
                i=i+1;
                labels = varargin{i};
        end
    end
end

% ensure field names have standard format
for i=1:numel(D)
   D(i).fields = standardize_fieldnames(D(i).fields); 
end

% trim day information from fieldnames
% strip day info from non-circadian fields
for i=1:numel(D)
    [a,m] = parse_fieldnames(D(i).fields);
    is_circ = strcmp(a,'Circadian');
    D(i).fields(~is_circ) = cellfun(@(a,m) sprintf('%s %s',a,m),a(~is_circ),m(~is_circ),'UniformOutput',false);
end

% get unique field names
field_list = unique(cat(1,D.fields));

% set plotting parameters
height = 0.3;
width = 0.5;
vspace = 0.3;
hspace = 0.105;


% iterate over fields
nRows = 13;
nCols = 10;
colors = {'b',[1 0.5 0],'m'};
ax_opts = {'XLim',[-4 4],'YLim',[0 0.8],'TickLength',[0 0],...
    'XTick',[],'YTick',[],'Units','inches','FontSize',6};
%%
for i=1:numel(field_list)
    
    % open new figure if necessary
    if mod(i-1,nRows*nCols)==0
       figure('Units','inches','Position',[.5 .5 7.5 9]);
       ax = subplot_array(nRows,nCols,ax_opts{:});
       drawnow;
    end
    plot_num = mod(i-1,nRows*nCols)+1;
    ah = ax(plot_num);
    hold(ah,'on');
    
    % get current field
    f = field_list{i};
    
    % iterate over data structs
    for j=1:numel(D)
        % find matching data field
        fidx = find(strcmpi(D(j).fields,f),1);
        fdat = D(j).data(:,fidx);
        if any(~isnan(fdat))
            plot_kde(fdat,ah,[-4 4],colors{j});
        end
        txt = lower(sprintf('n=%i',sum(~isnan(fdat))));
        if mod(j,2)
            th = text(.25,ah.YLim(2)*0.85,txt,'Color',colors{j},'FontSize',6,'Parent',ah);
        else
            th = text(.25,ah.YLim(2)*0.60,txt,'Color',colors{j},'FontSize',6,'Parent',ah);
        end
    end
    if mod(plot_num,nCols)==1
        ylabel(ah,'density');
        set(ah,'YTick',[0 .4 .8]);
    end
    if ceil(plot_num/nCols)==nRows
        set(ah,'XTick',[-4 0 4]);
    end
    [a,m,d] = parse_fieldnames({f},true);
    m{1}(m{1}=='_') = ' ';
    if ~isnan(d)
        x_label = {a{1};sprintf('%s (%i)',m{1},d)};
    else
        x_label = {a{1};m{1}};
    end
    xlabel(ah,x_label); 
    ah.Position(1) = (mod(plot_num-1,nCols))*(width+hspace) + 0.5;
    ah.Position(2) = (nRows - ceil(plot_num/nCols))*(height+vspace) + 0.5;
    ah.Position(3:4) = [width height];
    drawnow;
end
