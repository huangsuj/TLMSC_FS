function sigma = optSigma(X)
    N = size(X,1); %N���������ͼ��
    dist = EuDist(X,X);  %����ŷʽ�������
    dist = reshape(dist,1,N*N); %��Ƕ������reshape�����������µ��������������������ά��,����Ԫ�ظ�������
    sigma = median(dist); %ȡÿһ�еĴ�С�ŵ����м�ֵ���������ȡ��ֵ