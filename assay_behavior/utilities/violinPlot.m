function out=violinPlot(data, varargin)

bufferWidth=1;
numPts=100;
violinColor=[0.8 0.8 0.8];

switch class(data)
    case 'double'
        numVectors=size(data,2);
        dataTemp={};
        for i=1:numVectors
            dataTemp=[dataTemp data(:,i)]; 
        end
        data=dataTemp;
        means = mean(data);
    case 'cell'
        numVectors=length(data);
        [means,~,ci95s,~] = cellfun(@(x) normfit(x(~isnan(x))), data, ...
            'UniformOutput', false);
end

N = cellfun(@(d) sum(~isnan(d)), data);
if size(N,1) < size(N,2)
   N = N'; 
end
labels = arrayfun(@(i) sprintf('group #%i',i), 1:numel(N), 'UniformOutput', false);
for i=1:numel(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'Labels'
                i=i+1;
                labels = varargin{i};
        end
    end
end

labels = cellfun(@(s,n) sprintf('%s (n=%i)',s,n),...
    labels, num2cell(N), 'UniformOutput', false);

hold on;

fs=zeros(numPts,numVectors);
xis=zeros(numPts,numVectors);
fs2=zeros(numPts+2,numVectors);
xis2=zeros(numPts+2,numVectors);

for i=1:numVectors
    if any(~isnan(data{i}))
        [fs(:,i),xis(:,i)] = ksdensity(data{i});    
        xis2(:,i)=sort([xis(:,i);0;1]);
        fs2(:,i)=interp1(xis(:,i),fs(:,i),xis2(:,i),'spline',0);
    end
end

fs=fs/(max(max(fs))*(2+bufferWidth));


for i=1:numVectors
    handles.violins(i)=fill([i-fs(:,i)' i+fliplr(fs(:,i)')],[xis(:,i)' fliplr(xis(:,i)')],violinColor);
    handles.violins(i).EdgeColor='none';
end

vx = [repmat(1:numVectors,2,1); NaN(1,numVectors)];
vy = [cat(2,ci95s{:}); NaN(1,numVectors)];
handles.ci95s = plot(vx(:),vy(:),'Color',[.45 .45 .45],'LineWidth',1.5);
% handles.means = plot(1:numVectors,cat(2,means{:}),'Marker','o','MarkerFaceColor','k',...
%     'LineStyle','none','MarkerSize',4,'MarkerEdgeColor','none');
vx = [repmat(1:numVectors,2,1)+repmat([-.1;.1],1,numVectors); NaN(1,numVectors)];
vy = [repmat(cat(2,means{:}),2,1); NaN(1,numVectors)];
handles.means = plot(vx(:),vy(:),'Color','k','LineWidth',.5);

set(gca,'XTick',1:numel(data),'XTickLabels',labels,'XTickLabelRotation',90);

out.xis=xis;
out.fs=fs;
out.handles=handles;

