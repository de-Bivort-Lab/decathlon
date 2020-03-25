%% Get parent directory of all decathlon files
load('D:\Decathlon Raw Data\decathlon 1-2019\meta\culling_ID_permutation.mat');
fDir = 'D:\Decathlon Raw Data\decathlon 1-2019\data\';

    


%% create data matrix and create labels array

% query max number of data points
nFields = 0;
for i = 1:length(dec)
    dec(i).n = length(dec(i).ID);
    nFields = nFields + length(dec(i).fields);
end

cnFields=0;
for i = 1:length(circ)
    circ(i).n = length(circ(i).ID);
    cnFields = cnFields + length(circ(i).fields);
end

% model nuisance variable effects and replaced with residuals if necessary
%[dec,declm] = modelEffects(dec,nFields,'TimeofDay',false);
%[circ,circlm] = modelEffects(circ,cnFields,'TimeofDay',false);

% initialize data mat (observations x variables)
nFields = nFields + cnFields;
n={dec(:).ID};
n=cat(1,n{:});
dMat = NaN(max(n),nFields);
dFields = cell(nFields,1);
fct = 0;
nDays = max([dec.day circ.day]);

for i=1:nDays
    
    if i<=length(circ) && ~isempty(circ(i).fields)
        
        f = circ(i).fields;

        for j=1:length(f)

            fct=fct+1;
            dMat(circ(i).ID,fct) = circ(i).data.(f{j});
            dMat(circ(i).ID(~circ(i).data.filter),fct) = NaN;
            dFields(fct) = {[circ(i).name ' ' f{j} ' (' num2str(i) ')']};
            
            if any(dMat(:,end)==0)
                disp(i);
            end
        end 
    end
    
    if i<=length(dec) && ~isempty(dec(i).fields)
        
        f = dec(i).fields;

        for j=1:length(f)

            fct=fct+1;
            dMat(dec(i).ID,fct) = dec(i).data.(f{j});
            dMat(dec(i).ID(~dec(i).data.filter),fct) = NaN;
            dFields(fct) = {[dec(i).name ' ' f{j} ' (' num2str(i) ')']};
        end 
    end
end

% delete empty struct entries
dec(arrayfun(@(x) isempty(dec(x).data),1:length(dec)))=[];
circ(arrayfun(@(x) isempty(circ(x).data),1:length(circ)))=[];


dMat(~any(~isnan(dMat')),:)=[];
clearvars -except n dec circ dMat dFields nFields declm circlm

% initialize Decathlon data struct
D.data = dMat;
D.fields = dFields;

%% unfiltered raw data

% plot number of active flies on each day of testing
plot_numActive(D);

% plot correlation matrix and accompanying p-values
[fh,r,p]=plotCorr(D.data,'Labels',D.fields,'Cluster',true,...
    'Signed',false,'PlotTitle','Decathlon 3');

% plot pairwise sample size
plot_pairwiseSampleSize(D);

%% separate measures into distinct clusters with no apriori hypothesis of correlation

D_col = collapseMetrics(D,'Fields','all','Mode','PCA');

% plot correlation matrix for collapsed matrix and accompanying p-values
[fh,r,p]=plotCorr(D_col.data,'Labels',D_col.fields,'Cluster',true,...
    'Signed',true,'PlotTitle','Decathlon 3 - apriori collapsed');


%% plot correlations sorted by date

% collapse circadian metrics
D_col = collapseMetrics(D,'Fields','circadian','Mode','average');
plotCorr_byDate(D_col.data,D_col.fields);

%% plot correlation matrices per group

groups = [{{'circling';'right_bias'};{'speed';'nTrials'};{'clumpiness'};{'switchiness'};...
    {'Circadian'}}; unique(arrayfun(@(d) d.name, dec, 'UniformOutput', false));...
    arrayfun(@(i) sprintf('(%i)',i), unique(cat(1,circ.day)), 'UniformOutput', false)];

corrMat_byGroup(groups, dFields, fDir);


