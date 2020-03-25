function out=decathlonLasso(X,Y)

numPredictors=size(X,2);
numDataPoints=length(Y);

lambda=0.4;

CVfold=10;
[B,S]=lasso(X,Y,'CV',CVfold);
lambda=S.LambdaMinMSE;
out.lambda=lambda;
out.Bfull=B(:,S.IndexMinMSE);
%out.Sfull=S;


h=waitbar(0,'progress');
numResamples=100;
leaveOneOut=1;
if leaveOneOut==1
    numResamples=numDataPoints;
    out.preds=zeros(numDataPoints,1);
    out.actuals=zeros(numDataPoints,1);
end


out.rs=zeros(numResamples,1);
%out.dfs=zeros(numResamples,1);
out.Bs=NaN(size(X,2),numResamples);
%out.Bsconcatenated=[];

for i=1:numResamples
    
    if leaveOneOut==1
        trainPts=1:numDataPoints;
        trainPts(i)=[];
        [B,S]=lasso(X(trainPts,:),Y(trainPts),'Lambda',lambda);
        out.preds(i)=X(i,:)*B;
        out.actuals(i)=Y(i);
        out.Bs(:,i)=B;
        
    else
        testPts=randperm(numDataPoints,round(numDataPoints/CVfold));
        trainPts=1:numDataPoints;
        trainPts(testPts)=[];
        
        
        [B,S]=lasso(X(trainPts,:),Y(trainPts),'Lambda',lambda);
        
        out.Bs(:,i)=B;
        %out.Bsconcatenated=[out.Bsconcatenated;find(B~=0)];
        out.rs(i)=corr2(X(testPts,:)*B,Y(testPts));
        %out.dfs(i)=S.DF;
    end
    
    waitbar(i/numResamples,h);
end

[out.rs,out.ps] = corr(out.preds,out.actuals);

close(h);