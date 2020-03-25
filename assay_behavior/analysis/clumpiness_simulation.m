
div = [1000; 100; 10; 2];

for j=1:4
    nreps = 5;
    n = 500;
    L = 50000;
    clumpiness = NaN(nreps,n);
    bout_num = NaN(nreps,n);
    max_c = L/div(j);
    stds_iti = NaN(nreps,n);

    for i=1:nreps
        fprintf('iteration %i of %i\n',i,nreps);
       n_c = randi(max_c,[n 1]);
       rr = num2cell(rand(L,n),1);
       [~,rr_sort] = cellfun(@(r) sort(r), rr, 'UniformOutput', false);
       rr_sort = cellfun(@(r,nn) sort(r(1:nn*2)), rr_sort, num2cell(n_c)', 'UniformOutput', false);
       idx = cellfun(@(r) [r(1:2:end) r(2:2:end)], rr_sort, 'UniformOutput', false);
       nbouts = num2cell(cellfun(@(ii) size(ii,1), idx));
       iti = cellfun(@(ii) diff([[0;ii(1:end-1,2)] ii(:,1)],1,2), idx, 'UniformOutput', false);
       clumpiness(i,:) = cellfun(@(ii,nn) std(ii)/(L/(nn+1)), iti, nbouts);
       
       bout_num(i,:) = cat(2,nbouts{:});
    end

    % calculate mean at each number
    bout_num = bout_num(:);
    clumpiness = clumpiness(:);
    mean_line = NaN(max_c,1);
    for i=1:max_c
        mask = bout_num==i;
        mean_line(i) = mean(clumpiness(mask));
    end
    % plot results
    subplot(2,2,j);
    % PCARegressionCI(nanzscore([clumpiness(:),bout_num(:)]),...
    %     [false(numel(clumpiness),1) false(numel(clumpiness),1)],gca);
    pretty_scatter(clumpiness(:),bout_num(:),'k');
    hold on; plot(mean_line,1:max_c,'r-','LineWidth',1);
    xlabel('clumpiness');
    ylabel('bout number');
    rs = corr([clumpiness(:),bout_num(:)],'Type','Spearman');
    rp = corr([clumpiness(:),bout_num(:)],'Type','Pearson');
    title(sprintf('max choices = %i; spearman r=%0.2f; pearson r=%0.2f',max_c,rs(1,2),rp(1,2)));
end

%%

nreps = 25;
n = 500;
stds_iti = NaN(nreps,n);
mean_iti = NaN(nreps,n);
clumpiness = NaN(nreps,n);
bout_num = NaN(nreps,n);
xx = round(logspace(log10(1),log10(25000),nreps));


for i=1:nreps
   fprintf('iteration %i of %i\n',i,nreps);
   %n_c = randi(i,[n 1]);
   n_c = repmat(xx(i),[n 1]);
   rr = num2cell(rand(L,n),1);
   [~,rr_sort] = cellfun(@(r) sort(r), rr, 'UniformOutput', false);
   rr_sort = cellfun(@(r,nn) sort(r(1:nn*2)), rr_sort, num2cell(n_c)', 'UniformOutput', false);
   idx = cellfun(@(r) [r(1:2:end) r(2:2:end)], rr_sort, 'UniformOutput', false);
   nbouts = num2cell(cellfun(@(ii) size(ii,1), idx));
   iti = cellfun(@(ii) diff([[0;ii(1:end-1,2)] ii(:,1)],1,2), idx, 'UniformOutput', false);
   clumpiness(i,:) = cellfun(@(ii,nn) std(ii)/(L/(nn+1)), iti, nbouts);
   bout_num(i,:) = cat(2,nbouts{:});
   stds_iti(i,:) = cellfun(@(ii) std(ii), iti);
   mean_iti(i,:) = cellfun(@(nn) L/(nn+1), nbouts);
end

subplot(2,1,1);
lh1 = plot(xx,mean(stds_iti,2),'r-','Linewidth',1.5); 
hold on; 
lh2 = plot(xx,mean(mean_iti,2),'b-','Linewidth',1.5);
ylabel('iti');
xlabel('choice number');
legend([lh1 lh2],{'std(iti)','mean(iti)'});
subplot(2,1,2);
plot(xx,mean(stds_iti,2)./mean(mean_iti,2),'k-','Linewidth',1.5);
ylabel('std(iti)/mean(iti)');
xlabel('choice number');



for i=1:nreps
    ah = subplot(5,5,i);
        % calculate mean at each number
    mean_line = NaN(i,1);
    mu = NaN(i,1);
    for j=1:i
        mask = bout_num(i,:)==j;
        mean_line(j) = mean(stds_iti(i,mask));
        mu(j) = mean(mean_iti(i,mask));
    end
    plot(mean_line,mu,'r-','LineWidth',1);
    xlabel('std(iti)');
    ylabel('mean(iti)');
end
