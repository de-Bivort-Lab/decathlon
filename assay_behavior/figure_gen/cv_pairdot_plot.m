function cv_pairdot_plot(A,B,varargin)

labels = {'A';'B'};
for i=1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'Labels'
                i = i+1;
                labels = varargin{i};
        end
    end
end

% match and sort data fields across decathlon data sets
[A,B] = pairFields(A,B, 'Trim', true);

% compute coefficient of variance for each metric in each data set
A_cov = cellfun(@(d) get_cov(d), num2cell(A.data,1));
B_cov = cellfun(@(d) get_cov(d), num2cell(B.data,1));

% initialize plotting vertices for lines connecting paired dots and plot
ticks = [0.5;1.5];
vx = repmat([ticks;NaN],numel(A_cov),1);
vy = [A_cov; B_cov; NaN(size(A_cov))];
vy = vy(:);
plot(vx, vy, 'Color', [.3 .3 .3], 'LineWidth', .25);

% set plotting options for dots and plot
opts = {'Marker'; 'o'; 'LineStyle'; 'none';...
        'MarkerFaceColor'; 'k'; 'MarkerEdgeColor'; 'none';...
        'MarkerSize'; 4; 'LineWidth'; 1};
hold on;
plot(ones(numel(A_cov),1).*ticks(1), A_cov, opts{:});
plot(ones(numel(B_cov),1).*ticks(2), B_cov, opts{:});
plot(ticks,[mean(A_cov),mean(B_cov)],'r','LineWidth',1);
ph = plot(ticks,[mean(A_cov),mean(B_cov)],opts{:});
ylabel('Coefficient of Variation');
xlabel('Decathlon Group');
set(ph,'MarkerSize',6,'MarkerFaceColor','r');
set(gca,'XLim',[0 2],'XTick',ticks,'XTickLabel',labels);


% compute coefficient of variation for data vector
function cov = get_cov(data)
data = data(~isnan(data));
cov = std(data)/mean(abs(data));