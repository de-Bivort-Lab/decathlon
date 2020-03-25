
fdir = ['D:\decathlon_preprint_code_data_figures\decathlon_analysis\'...
    'matrices\decathlon_paper\decathlon_final\'];
D = load_decathlon_structs(fdir,'D13_als');
%D = pair_decathlon_structs(D,'CollapseMode','PCA','CollapseFields','all');
d_label = {'inbred';'outbred'};

figure; 
hold on;
lhs = gobjects(2,1);
colors = {[0 0 0];[1 0 0]};


nreps = 1000;
bin_cts = cell(2,1);
bin_cts(:) = {NaN(nreps)};

for i=1:2
    d = D(i).data;
    
    % plot connected components output
    for qq=1:nreps
        idx = randi(size(d,1),[size(d,1) 1]);
        if ~mod(qq,20)
            fprintf('bootstrap iter %i of %i\n',qq,nreps);
        end
        [~,bin_cts{i}(qq,:)]=decathlonConnCompSweep(d(idx,:),nreps);
    end
    x = 0:20:size(d,2);
    cts = cellfun(@(c) histc(c,x),...
        num2cell(bin_cts{i},2),'UniformOutput',false);
    cts = cat(1,cts{:});
    mu = mean(cts);
    lb = prctile(cts,2.5);
    ub = prctile(cts,97.5);
    vx = [1 x fliplr(x)];
    vy = [ub(1) lb fliplr(ub)];
    vy = log10(vy(:));
    vy(vy<-3)=-3;
            
    lhs(i) = plot(x,log10(mu),'Color',colors{i},'Linewidth',1); hold on;
    patch('XData',vx(:),'YData',vy(:),'FaceColor',colors{i},'FaceAlpha',0.5);
    set(gca,'YLim',[0 3],'XLim',[1 x(end)]);
    ylabel('log(count)');
    xlabel('effective dimensionality');
end
legend(lhs,d_label);
title('Decathlon Behavior');