%% Get parent directory of all decathlon files

fDir = uigetdir('D:\decathlon_data_and_analysis','Select parent directory containing experiment data');
save_dir = uigetdir('D:\decathlon_data_and_analysis','Select directory to save output figures');
exp_name = 'D2';
ID = [];
%load('D:\decathlon_data_and_analysis\decathlon 1-2019\meta\culling_ID_permutation.mat');
options = {'FieldFunction',@getDataFields_legacy,'CullingPermutation',ID};

%% build the decathlon data struction from individual assay files

D = build_DecathlonDataStruct(fDir,options{:});

%% unfiltered raw data

% plot number of active flies on each day of testing
fh = plot_numActive(D);
savefig(fh,sprintf('%s/%s%s',save_dir,exp_name,'_numactive'));

% plot correlation matrix and accompanying p-values
[fh,r,p]=plotCorr(D.data,'Labels',D.fields,'Cluster',false,...
    'Signed',false,'PlotTitle',sprintf('%s - all metrics',exp_name),...
    'SavePath',sprintf('%s/%s%s',save_dir,exp_name,'_allmetrics'));

% plot pairwise sample size
[fh, D.pw_sample_size] = plot_pairwiseSampleSize(D);
savefig(fh,sprintf('%s/%s%s',save_dir,exp_name,'_pairwise_sampling'));

%% separate measures into distinct clusters with no apriori hypothesis of correlation

D_col = collapseMetrics(D,'Fields','all','Mode','PCA');

% plot correlation matrix for collapsed matrix and accompanying p-values
[fh,r,p]=plotCorr(D_col.data,'Labels',D_col.fields,'Cluster',false,...
    'Signed',true,'PlotTitle',sprintf('%s - apriori collapsed',exp_name),...
    'SavePath',sprintf('%s/%s%s',save_dir,exp_name,'_apriori'));


%% plot correlations sorted by date

% collapse circadian metrics
D_col = collapseMetrics(D,'Fields','circadian','Mode','average');
fh = plotCorr_byDate(D_col.data,D_col.fields,...
    'PlotTitle',sprintf('%s - assays by day',exp_name),...
    'SavePath',sprintf('%s/%s%s',save_dir,exp_name,'_assaysbyday'));

%% plot correlation matrices per group

% define plotting groups
assay_names = regexp(D.fields,'\<[A-Z][\w|-]*','match');
assay_names = cellfun(@(s) s{1}, assay_names, 'UniformOutput', false);
days = regexp(D.fields,'(?<=\()[0-9]*(?=\))','match');
days = cellfun(@(s) s{1}, days, 'UniformOutput', false);
groups = [{{'circling';'right_bias'};{'speed';'nTrials'};{'clumpiness'};...
    {'switchiness'};{'Circadian'}}; {unique(assay_names)}; {unique(days)}];

% plot by group
corrMat_byGroup(groups, D, fDir);


