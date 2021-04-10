%% simulate several conn comp plots

% define data dimensions
m = 400;        % number of observations
n = 30;        % number of features
q = [2 6];         % ground truth dimensionality
p = [3 5];          % number of features in each correlated cluster
cluster_cov = .4;
fhs = gobjects(3,1);
for i=1:numel(fhs)
   fhs(i) = figure; 
end

for ii=1:numel(q)
    for jj=1:numel(p)
            
        ct = (ii-1)*numel(q) + jj;

        % impart correlations to cluster by adding shared noise
        sim_cov = zeros(n);
        for i=1:q(ii)
            j = (i-1)*p(jj);
            if i > q(ii)/2
                sim_cov(j+1:j+p(jj),j+1:j+p(jj)) = cluster_cov/2;
            else
                sim_cov(j+1:j+p(jj),j+1:j+p(jj)) = cluster_cov;
            end
        end
        sim_cov(logical(diag(ones(n,1))))=1;
        sim_data = mvnrnd(zeros(n,1),sim_cov,m);
        sim_data = zscore(sim_data);

        % display covariance matrix
        figure(fhs(1));
        ah=subplot(2,2,ct);
        imagesc(corr(sim_data));
        k = n-q(ii)*(p(jj)-1);
        title({sprintf('num features = %i, clusters = %i',n,q(ii));...
            sprintf('effective dimensions = %i',k)});
        axis('equal','tight');
        colormap(ah,nanticoke);
        colorbar;
        caxis([-1 1]);
        ah.Units = 'inches';
        ah.Position(3:4) = .85;
        drawnow;
        

        % plot connected components output
        out = decathlonConnCompDropN(sim_data);
        figure(fhs(2));
        ah=subplot(2,2,ct);
        imagesc(out);
        hold on;
        plot([k k],[0.5 size(out,1)+.5],'k--');
        ylabel('features dropped');
        xlabel('effective dimensionality');
        title('connected components heatmap');
        axis('equal','tight');
        colormap(ah,logjet_cmap);
        colorbar; caxis([0 1]);
        ah.Units = 'inches';
        ah.Position(3:4) = .85;
        drawnow;

        % plot connected components output
        nreps = 100;
        bin_cts = zeros(1,n);
        for qq=1:nreps
            idx = randi(size(sim_data,1),[m 1]);
            if ~mod(qq,20)
                fprintf('bootstrap iter %i of %i\n',qq,nreps);
            end
            bin_cts = bin_cts + ...
                histc(conn_comp_spectrum(corr(sim_data(idx,:)),nreps),1:n);
        end

        figure(fhs(3));
        ah=subplot(2,2,ct); hold on;
        bar(1:n,log10(bin_cts),'FaceColor',[.65 .65 .65],'EdgeColor','none');
        peaks = [(n-q(ii)*p(jj))+q(ii), (n-q(ii)/2*p(jj))+q(ii)/2, n];
        for pp = 1:numel(peaks)
           plot([peaks(pp) peaks(pp)],[0 4],'k--');
        end
        ylabel('log(count)');
        xlabel('connected components');
        ah.Units = 'inches';
        ah.Position(3:4) = .85;
        drawnow;
    end
end


%% Plot Scree plot derivative

% define data dimensions
m = 5000;        % number of observations
n = 30;        % number of features
q = [3 6];         % ground truth dimensionality
p = [2 5];          % number of features in each correlated cluster
cluster_cov = [0.25 0.5];
fhs = gobjects(2,1);
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
            %sim_data = mvnrnd(zeros(n,1),sim_cov,m);
            sim_data = mvnrnd(zeros(size(cov_mat,1),1),cov_mat,m);
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
            
            figure(fhs(2));
            [~,~,~,~,var_exp] = pca(sim_data);
            ah=subplot(4,2,ct); hold on;
            bar(-diff(var_exp),'FaceColor',[.75 .75 .75]);
            plot([k k],[0 5],'k--');
            ylabel({'-\Delta Variance';'Explained (%)'});
            xlabel('PC. rank');
            set(ah,'YLim',[0 5],'TickLength',[0 0]);
          
        end
    end
end

%% NEW FIGS9 analyses


% 1. repeat figs9 - 1) shuffle data and produce new cov matrix 
%    2) take empirical distribution of correlations and substract 
%    from it the estimate of the noise distribution 3) populate matrix 
%    with noise matched to shuffle noise, then overlay on top of values 
%    generated from (2)
% 2. Take out "correlated clusters" and take out "cluster covariance"
%    have 3 clusters at lower level and 3 at the higher covariance
% 3. Rows to keep - 1,4,6,7,8
% 4. Change num clusters to 2 and 6 (rather than 3 and 6)


D = load_decathlon_structs(pwd,'D_als_filled_batch_merged');
D_p = pair_decathlon_structs(D,'ImputeMode','none','CollapseFields','none');

%%
bins = linspace(-1,1,1000);
nfeat = size(D_p(1).data,2);
fh = figure;
nreps = 100;
nflies = size(D_p(i).data,1);
for i=1:numel(D_p)
    % compute correlation matrix
    r_obs = corr(D_p(i).data);
    uidx = upper_triangle_idx(length(r_obs));
    r_obs = r_obs(uidx);
    r_shuf = corr(shuffle_columns(D_p(i).data));
    r_shuf = r_shuf(uidx);

    % compute kde
    obs_dist = ksdensity(r_obs,bins);
    shuf_dist = ksdensity(r_shuf,bins);
    diff_dist = obs_dist - shuf_dist;
    diff_dist(diff_dist<0) = 0;
    
    % generate a corr mat from random distribution
    shuf_cdf = cumsum(shuf_dist)./sum(shuf_dist);
    [row,col] = ind2sub([nfeat nfeat],uidx);
    
    % generate correlation values from signal distribution
    nsignal_feat = floor((sum(diff_dist)/sum(obs_dist))*numel(uidx));
    diff_cdf = cumsum(diff_dist)./sum(diff_dist);


    comp_cts = zeros(nreps);
    for j=1:nreps
        if mod(j,10)==0
           fprintf('Conn. comp. simulation %i of %i\n',j,nreps); 
        end
        % generate simulated corr mat noise
        rand_draws = rand(size(uidx));
        [~,min_idx] = min(abs(shuf_cdf - rand_draws),[],2);
        r_sim = diag(ones(nfeat,1));
        r_sim(uidx) = bins(min_idx);
        r_sim(sub2ind([nfeat nfeat],col,row)) = bins(min_idx);
        
        % add signal
        rand_draws = rand([nsignal_feat 1]);
        [~,min_idx] = min(abs(diff_cdf - rand_draws),[],2);
        rand_uidx = uidx(randperm(numel(uidx),nsignal_feat));
        [sig_row,sig_col] = ind2sub([nfeat nfeat],rand_uidx);
        r_sim(rand_uidx) = bins(min_idx);
        r_sim(sub2ind([nfeat nfeat],sig_col,sig_row)) = bins(min_idx);
        
        % compute conn comp spectrum
        comp_cts(j,:) = conn_comp_spectrum(r_sim,nreps);
    end
   
    % plot conn comp spetrogram
    ah = subplot(3,2,i+2);
    cts = histc(comp_cts(:),1:nfeat);
    bar(1:nfeat,log10(cts),'FaceColor',[.65 .65 .65],'EdgeColor','none');
    xlabel('conn. comp');
    ylabel('log(counts)');
    if mod(i,2)
        title('Inbred seeded');
    else
        title('Outbred seeded');
    end
    drawnow;
    ah.Units = 'inches';
    ah.Position(3:4) = [2 1];
    drawnow;
    
    
    % plot matrix
    fh = figure;
    z = linkage(r_sim);
    [~,~,p]=dendrogram(z,0);
    close(fh);
    ah = subplot(3,2,i);
    imagesc(r_sim(p,p));
    colormap(nanticoke);
    caxis([-1 1]);
    colorbar;
    set(ah,'XTick',[],'YTick',[]);
    if mod(i,2)
        title('Inbred seeded');
    else
        title('Outbred seeded');
    end
    drawnow;
    ah.Units = 'inches';
    axis('equal','tight');
    ah.Position(3:4) = 1;
    drawnow;
    
    comp_cts = zeros(nreps);
    for j=1:nreps
        if mod(j,10)==0
           fprintf('Conn. comp. simulation %i of %i\n',j,nreps); 
        end
        d = D_p(i).data(randi(nflies,[nflies 1]),:);
        
        % compute conn comp spectrum
        comp_cts(j,:) = conn_comp_spectrum(corr(d),nreps);
    end
    
    ah = subplot(3,2,i+4);
    cts = histc(comp_cts(:),1:nfeat);
    bar(1:nfeat,log10(cts),'FaceColor',[.65 .65 .65],'EdgeColor','none');
    xlabel('conn. comp');
    ylabel('log(counts)');
    if mod(i,2)
        title('Inbred');
    else
        title('Outbred');
    end
        drawnow;
    ah.Units = 'inches';
    ah.Position(3:4) = [2 1];
    drawnow;
end
