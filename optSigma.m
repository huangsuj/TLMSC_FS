function sigma = optSigma(X)
    N = size(X,1); %N是输入的视图数
    dist = EuDist(X,X);  %计算欧式距离矩阵
    dist = reshape(dist,1,N*N); %内嵌函数：reshape函数用于重新调整矩阵的行数、列数、维数,但是元素个数不变
    sigma = median(dist); %取每一列的从小排到大中间值，再排序后取中值