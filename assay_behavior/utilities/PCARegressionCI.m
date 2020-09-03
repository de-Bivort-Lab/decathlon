function out=PCARegressionCI(data,ah,varargin)
% performs linear regression on the two column variables of data
% using PCA to find the regression fit, i.e., using normals to compute
% residuals
%
% performs bootstrap resampling of this regression to determine the
% confidence interval of the regression line.

r=3.5;
xlim = [-r r];
ylim = [-r r];
plotBool=true;     % make a figure at the end?
for i=1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'XLim'
                i=i+1;
                xlim = varargin{i};
            case 'YLim'
                i=i+1;
                ylim = varargin{i};
            case 'Plot'
                i=i+1;
                plotBool = varargin{i};
        end
    end
end

% configure parameters here:
numReps=1000;   % # of bootstrap replicates
CI=95;          % confidence interval to be estimated


numPts=size(data,1);

% generate points for line and CI
numXs=3;
minX=min(data(:,1));
maxX=max(data(:,1));

xPts=linspace(xlim(1),xlim(2),numXs);

out.fitVals=zeros(numReps,numXs);
obs = pca(data);
x_bar=mean(data(:,1));
y_bar=mean(data(:,2));
m=obs(2,1)/obs(1,1);
fx = @(x) m*x + (y_bar - m*x_bar);
obs =m*xPts+(y_bar-m*x_bar);

% bootstrap loop
m_bs = NaN(numReps,1);
for i=1:numReps
    dataTemp=data(randperm(numPts,round(numPts*0.7)),:);
    [eigenvectors,~,~] = pca(dataTemp);
    x_bar=mean(dataTemp(:,1));
    y_bar=mean(dataTemp(:,2));
    m=eigenvectors(2,1)/eigenvectors(1,1);
    out.fitVals(i,:)=m*xPts+(y_bar-m*x_bar);
    m_bs(i) = m;
end

%out.fitVals = nanzscore(out.fitVals);
m = median(m_bs);
fx = @(x) m*x + (y_bar - m*x_bar);
obs =m*xPts+(y_bar-m*x_bar);

% save CI bounds
out.lower=prctile(out.fitVals,(100-CI)/2);
out.upper=prctile(out.fitVals,(100-(100-CI)/2));

y = fx(xPts([1 end]));
y = round(y*100)/100;
if(any(y > 3.5 | y < -3.5))
    x = linspace(-15,15,10000);
    [~,x1] = min(abs(fx(x)+3.5));
    [~,x2] = min(abs(fx(x)-3.5));
    out.fit_x = sort(x([x1 x2]));
    out.fit_y = fx(out.fit_x);
else
    out.fit_y = y;
    out.fit_x = xPts([1 end]);
end



% make a figure if plotBool==1.
if plotBool==1
    hold on;
    opts = {'Marker'; 'o'; 'LineStyle'; 'none';...
    'MarkerEdgeColor'; 'none';'MarkerSize'; 1.5; 'LineWidth'; 1};
    patch([xPts fliplr(xPts)],[out.lower fliplr(out.upper)],[0 0 0],...
        'Parent',ah,'faceAlpha',0.2,'edgeColor','none');
%     plot([xPts(1) xPts(end)],[median(out.fitVals(:,1)) median(out.fitVals(:,end))],'k',...
%         'Parent',ah);
    plot([xPts(1) xPts(end)],[obs(1) obs(end)],'k','Parent',ah);
    plot(data(:,1),data(:,2),'Parent',ah,opts{:},...
        'MarkerFaceColor','k');
    set(ah,'XLim',xlim,'YLim',ylim);
end

    
