clear;clc;
addpath('tSVD','proxFunctions','solvers','twist','code_coregspectral');
addpath('ClusteringMeasure', 'LRR', 'Nuclear_norm_l21_Algorithm', 'unlocbox');
addpath('.\datasets');
dataset='NUS-WIDE'; 
load(dataset);
% truth = Y;     
numClus=length(unique(truth));
views=length(X);
dataNum=size(X{1},2);
for iv=1:views
    X{iv} = double(X{iv}');
end

MAXiter = 1000; 
REPlic = 20; % Number of replications for KMeans
repnum=10;                                     
M= numClus;  
alpha= 0.05; 
gamma = 0.001;
lambda = 0.0001;
beta = 0.2;
for i=1:repnum
        t1=clock;
        [A]=TLMSC_FS(alpha,gamma,lambda,beta,X,numClus);
        label=litekmeans(A, numClus, 'MaxIter', 100,'Replicates',REPlic);
        truth = double(truth);
        result=ClusteringMeasure(truth,label);

        t2=clock;
        tempACC(i) = result(1);
        tempNMI(i) = result(2);
        tempPurity(i) = result(3);
        tempARI(i) = result(4);
        tempFscore(i) = result(5);
        tempPrecision(i) = result(6);
        tempRecall(i) = result(7);
        Time(i) = etime(t2,t1);
        
        ACC = [mean(tempACC),std(tempACC)];
        ARI = [mean(tempARI),std(tempARI)];
        NMI = [mean(tempNMI),std(tempNMI)];
        Purity = [mean(tempPurity),std(tempPurity)];
        Fscore = [mean(tempFscore),std(tempFscore)];
        Precision = [mean(tempPrecision),std(tempPrecision)];
        Recall = [mean(tempRecall),std(tempRecall)];
        Time = [mean(Time),std(Time)]; 
 end

    ACC = roundn(ACC,-3)
    NMI = roundn(NMI,-3)
    Purity = roundn(Purity,-3)
    ARI = roundn(ARI,-3)
    Fscore =roundn(Fscore,-3)
    Precision =roundn(Precision,-3)
    Recall =roundn(Recall,-3)

    Time = roundn(Time,-3)
    





                

    
%%