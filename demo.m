clc;clear;
addpath C:\Users\24952\Desktop\CRAA
path = 'C:\Users\24952\Desktop\CRAA';
datapath = 'C:\Users\24952\Desktop\CRAA\';
datasetname = {'MSRCV1'};
for di = 1:length(datasetname)
    dataName = datasetname{di};
    load([datapath,'/',dataName,'.mat']);
    disp([' Current dataset  : ',dataName]);
    
    num_views = size(X,2);
    numClust = size(unique(gt),1);
    
    opts.V = size(X,2);
    opts.N = size(X{1},2);
    
    for i=1:opts.V
        X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1);
    end
    
    opts.SD = 0;
    M = [];
    
    for i=1:opts.V
        opts.SD = opts.SD + size(X{i},1);
        M = [M;X{i}];
    end
    
    Xopts.lambda_1=[0.001 0.01 0.1 1 10 100 1000 10000];
    Xopts.K = int_mod(opts.SD,numClust);
    
    idx = 1;
    idxq = 1;
    
    for li=1:length(Xopts.K)
        K=Xopts.K(li);
        for j=1:length(Xopts.lambda_1)
            lambda_1 = Xopts.lambda_1(j);
            for i=1:20
                [C,max_iter]=CRAA(X,opts,M,lambda_1,K);
                [ress(max_iter,:)] =  ClusteringMeasure(abs(C)+abs(C'),gt);
                res_comvsc = ress(end,:);
                disp(['K = ',num2str(K),'  lambda_1 = ',num2str(lambda_1), '  ACC = ', num2str(res_comvsc(:,2))]);
                idx = idx+1;
            end
        end
    end
end