function plot_tsne_position_samples(embedding)

firstFrame = randi(round(numel(embedding.z_speed{1})*0.8),[1 1]);
lastFrame = firstFrame + 100*20;
frame_range = firstFrame:lastFrame-4;
N = 6;
M = 8;
L = numel(embedding.z_data);
for i = 1:L
    
    if mod(i-1,M*N)==0
        figure;
    end
    
    subplot(M,N,mod(i-1,M*N)+1)
    hold on
    plot(embedding.z_data{i}(frame_range,1),'r','Linewidth',1);
    plot(embedding.z_data{i}(frame_range,2),'b','Linewidth',1);
    hold off
    set(gca,'YLim',[-100 100],'YTick',[],'XLim',[1 length(frame_range)],'XTick',[]);
    if mod(i,N)==1
        ylabel('position');
        set(gca,'YTick',[-100 0 100]);
    end
    if mod(i-1,M*N)+1>N*(M-1)
        set(gca,'Xtick',linspace(1,length(frame_range),5),'XTickLabel',0:5:20);
        xlabel('time(s)');
    end
    title(['fly #' num2str(i)],'FontWeight','normal','FontSize',6);
    drawnow;
end