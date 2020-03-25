function plot_mode_occupancy_pdf(dur)

all_dur = cat(1,dur{:});
all_dur = cat(1,all_dur{:});
all_dur = all_dur./100;
time_bins = logspace(log10(.01),log10(100),25);
dur_cts = histc(all_dur,time_bins);
dur_pdf = dur_cts./sum(dur_cts);
figure;
dur_idx = 1:numel(dur_cts);
dur_idx(dur_pdf==0)=[];
dur_pdf(dur_pdf==0)=[];
plot(dur_idx,log10(dur_pdf),'k-o','LineWidth',2);
set(gca,'XTick',1:6:numel(dur_cts),'XTickLabel',logspace(log10(0.01),log10(100),5),...
    'YTick',-8:0,'YTickLabel',logspace(log10(10^-8),0,9));
xlabel('occupancy duration (s)');
ylabel('PDF');