function [r,p] = corr_of_corrcoef(A,B, varargin)
% Compute and plot the correlation of r-values between decathlon data structures A and B

% set default settings and parse inputs
doplot = false;
titlestr='';

for i=1:numel(varargin)
	if ischar(varargin{i})
		switch varargin{i}
			case 'Plot'
				i = i + 1;
				doplot = varargin{i};
            case 'Title'
                i=i+1;
                titlestr = varargin{i};
		end
	end
end

% compute correlation matrix
[A.r,A.p] = corrcoef(A.data,'rows','pairwise');
[B.r,B.p] = corrcoef(B.data,'rows','pairwise');

% get the indices of the lower triangle of the correlation 
% matrix to avoid duplicate r-values and self-correlation
% r-values on the diagonal 
L=1:length(A.r);
subset = arrayfun(@(x) [L(L<x)' repmat(x,sum(L<x),1)],L,'UniformOutput',false);
subset = cat(1,subset{:});
subset = sub2ind(size(A.r),subset(:,1),subset(:,2));

% calculate correlation of r-values
[r,p] = corrcoef([A.r(subset),B.r(subset)],'rows','pairwise');

if doplot
    ah=gca;
    opts = {'Marker'; 'o'; 'LineStyle'; 'none';...
        'MarkerFaceColor'; 'k'; 'MarkerEdgeColor'; 'none';...
        'MarkerSize'; 2; 'LineWidth'; 2};
    plot(ah,A.r(subset),B.r(subset),'o',opts{:});
    disp(['r=' num2str(r(1,2)) ' (p=' num2str(p(1,2)) ')']);
    xlabel('r-value (dataset 1)');
    ylabel('r-value (dataset 2)');
    set(ah,'XLim',[-1 1],'YLim',[-1 1],'XTick',-1:0.5:1,'YTick',-1:0.5:1);
    text(ah,ah.XLim(1)+diff(ah.XLim)*0.075,ah.YLim(2)-diff(ah.YLim)*0.075,...
        ['r=' num2str(r(1,2),2) ' (p=' num2str(p(1,2),2) ')']);
    title(titlestr);
end
