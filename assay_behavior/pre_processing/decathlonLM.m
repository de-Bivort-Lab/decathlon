function out=decathlonLM(X,Y)

warning('off','stats:LinearModel:RankDefDesignMat');
numDataPoints=length(Y);

numResamples=1000;
out.rs=zeros(numResamples,1);
out.Bs=cell(numResamples,1);


CVfold=10;
PCsToUse=1:size(X,2);

wb=0;

if wb==1; h=waitbar(0); end

leave_one_out = true;
if leave_one_out
    
    ypred = NaN(numDataPoints,1);
    for i=1:numDataPoints

        testPts=i;
        trainPts=1:numDataPoints;
        trainPts(testPts)=[];

        rnaModel=fitlm(X(trainPts,PCsToUse),Y(trainPts));

        ypred(i) = predict(rnaModel,X(testPts,PCsToUse));

        out.Bs{i} = rnaModel.Coefficients{2:end,1};

        if wb==1;  waitbar(i/numResamples,h); end
    end
    
    out.rs=corr2(ypred,Y);
else
    for i=1:numResamples

        testPts=randperm(numDataPoints,round(numDataPoints/CVfold));
        trainPts=1:numDataPoints;
        trainPts(testPts)=[];

        rnaModel=fitlm(X(trainPts,PCsToUse),Y(trainPts));

        ypred = predict(rnaModel,X(testPts,PCsToUse));

        out.rs(i)=corr2(ypred,Y(testPts));
        out.Bs{i} = rnaModel.Coefficients{2:end,1};

        if wb==1;  waitbar(i/numResamples,h); end
    end
end

if wb==1;  close(h); end