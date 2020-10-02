%% simulate several conn comp plots

% define data dimensions
m = 500;        % number of observations
n = 30;        % number of features
q = [3 6];         % ground truth dimensionality
p = [2 5];          % number of features in each correlated cluster
cluster_cov = [0.25 0.5];
fhs = gobjects(3,1);
for i=1:numel(fhs)
   fhs(i) = figure; 
end

for ii=1:numel(q)
    for jj=1:numel(p)
        for kk=1:numel(cluster_cov)
            
            ct = (ii-1)*numel(q)*numel(p) + (jj-1)*numel(p) + kk;
            
            % impart correlations to cluster by adding shared noise
            sim_cov = zeros(n);
            for i=1:q(ii)
                j = (i-1)*p(jj);
                if i==q(ii)
                    sim_cov(j+1:j+p(jj),j+1:j+p(jj)) = cluster_cov(kk)/2;
                else
                    sim_cov(j+1:j+p(jj),j+1:j+p(jj)) = cluster_cov(kk);
                end
            end
            sim_cov(logical(diag(ones(n,1))))=1;
            sim_data = mvnrnd(zeros(n,1),sim_cov,m);
            sim_data = zscore(sim_data);

            % display covariance matrix
            figure(fhs(1));
            ah=subplot(4,2,ct);
            imagesc(cov(sim_data));
            k = n-q(ii)*(p(jj)-1);
            title({sprintf('num features = %i, clusters = %i',n,q(ii));...
                sprintf('effective dimensions = %i',k)});
            axis('equal','tight');
            colormap(ah,egoalley);
            colorbar;
            caxis([-1 1]);
            drawnow;

            % plot connected components output
            out = decathlonConnCompDropN(sim_data);
            figure(fhs(2));
            ah=subplot(4,2,ct);
            imagesc(out);
            hold on;
            plot([k k],[0.5 size(out,1)+.5],'k--');
            ylabel('features dropped');
            xlabel('effective dimensionality');
            title('connected components heatmap');
            axis('equal','tight');
            colormap(ah,logjet_cmap);
            colorbar; caxis([0 1]);
            drawnow;
            
            % plot connected components output
            nreps = 200;
            bin_cts = NaN(nreps);
            for qq=1:nreps
                idx = randi(size(sim_data,1),[size(sim_data,1) 1]);
                if ~mod(qq,20)
                    fprintf('bootstrap iter %i of %i\n',qq,nreps);
                end
                [~,bin_cts(qq,:)]=decathlonConnCompSweep(sim_data(idx,:),200);
            end
            x = 1:size(sim_data,2);
            cts = cellfun(@(c) histc(c,x)./numel(c),...
                num2cell(bin_cts,2),'UniformOutput',false);
            cts = cat(1,cts{:});
            mu = mean(cts);
            lb = prctile(cts,2.5);
            ub = prctile(cts,97.5);
            vx = [1 x fliplr(x)];
            vy = [ub(1) lb fliplr(ub)];
            
            figure(fhs(3));
            ah=subplot(4,2,ct);
            plot(x,mu,'k-','Linewidth',1); hold on;
            patch('XData',vx(:),'YData',vy(:),'FaceColor',[.5 .5 .5],'FaceAlpha',0.5);
            %histogram(bin_cts,1:size(sim_data,2),'Normalization','probability');
            
            hold on;
            plot([k k],[0 1],'k--');
            set(ah,'YLim',[0 1],'YTick',0:0.2:1);
            ylabel('probability');
            xlabel('effective dimensionality');
            drawnow;
        end
    end
end
