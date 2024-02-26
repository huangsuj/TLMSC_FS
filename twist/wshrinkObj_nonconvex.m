function [x,objV] = wshrinkObj_nonconvex(x,rho,sX, isWeight,mode)
e=10^-7;
if isWeight == 1
%     C = 2*sqrt(2)*sqrt(sX(3)*sX(2));
    C = sqrt(sX(3)*sX(2));
end
if ~exist('mode','var')
    % mode = 1是采用lateral slice的方法
    % mode = 2是采用front slice的方法
    % mode = 3是采用top slice的方法
    mode = 1;
end

X=reshape(x,sX);  %% sX=[n,n,k]
if mode == 1
    Y=X2Yi(X,3);
elseif mode == 3
    Y=shiftdim(X, 1);  %%将X维度往左循环移动1位，Y in [n,k,n]
else
    Y = X;
end
% 

Yhat = fft(Y,[],3);  %% 对Y进行傅里叶变换
% weight = C./sTrueV+eps;
% weight = 1;
% tau = rho*weight;
objV = 0;
if mode == 1
    n3 = sX(2);
elseif mode == 3   %% 张量维度变换之后，n3变为sX(1)
    n3 = sX(1);
else
    n3 = sX(3);
end

if isinteger(n3/2) 
    %% 偶数的时候执行以下语句
    endValue = int16(n3/2+1);
    for i = 1:endValue
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');  %%逐层进行SVD分解
        
        if isWeight   %% isWeight 取值为0，可以忽略
            weight = C./(diag(shat) + eps);
            tau = rho*weight;
            shat = soft(shat,diag(tau)); 
        else
            tau = rho;  %% 传的参数为1/rho
            % shat = max(shat - tau,0);   %% 由于Yhat已经进行了傅里叶变换，此处shat已进行了傅里叶变换
            % non-convex tensor rank minimization
            shat_Mf=diag(shat);
            shat_Gf=zeros(length(shat_Mf));
            for k=1:length(shat_Mf)
                triangle=(e-shat_Mf(k))^2-4*(tau-shat_Mf(k)*e);
                if triangle>=0
                    shat_Gf(k)=0.5*(shat_Mf(k)-e+sqrt(triangle));
                else  
                    shat_Gf(k)=0;
                end
            end 
            
        end                 
        shat=diag(shat_Gf);
        objV = objV + sum(shat(:));
        Yhat(:,:,i) = uhat*shat*vhat';
        if i > 1
            Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';  
            objV = objV + sum(shat(:));
        end
    end
    [uhat,shat,vhat] = svd(full(Yhat(:,:,endValue+1)),'econ');
    if isWeight
       weight = C./(diag(shat) + eps);
       tau = rho*weight;
       shat = soft(shat,diag(tau));
    else
       tau = rho;
       % shat = max(shat - tau,0);
       % non-convex tensor rank minimization
       shat_Mf=diag(shat);
       shat_Gf=zeros(length(shat_Mf));
       for k=1:length(shat_Mf)
           triangle=(e-shat_Mf(k))^2-4*(tau-shat_Mf(k)*e);
           if triangle>=0
               shat_Gf(k)=0.5*(shat_Mf(k)-e+sqrt(triangle));
           else
               shat_Gf(k)=0;
           end
       end 
       
    end
    shat=diag(shat_Gf);
    objV = objV + sum(shat(:));
    Yhat(:,:,endValue+1) = uhat*shat*vhat';
    
    %%奇数的时候执行以下语句
else
    endValue = int16(n3/2+1);  %int16(n3/2+1)
    for i = 1:endValue
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');
        if isWeight
            weight = C./(diag(shat) + eps);
            tau = rho*weight;
            shat = soft(shat,diag(tau));
            
        else
            tau = rho;
            %shat = max(shat - tau,0);  %%更新shat
            % non-convex tensor rank minimization
            shat_Mf=diag(shat);
            shat_Gf=zeros(length(shat_Mf));
            for k=1:length(shat_Mf)
                triangle=(e-shat_Mf(k))^2-4*(tau-shat_Mf(k)*e);
                if triangle>=0
                    shat_Gf(k)=0.5*(shat_Mf(k)-e+sqrt(triangle));
                else
                    shat_Gf(k)=0;
                end
            end
        end
        shat=diag(shat_Gf(k));
        objV = objV + sum(shat(:));
        Yhat(:,:,i) = uhat*shat*vhat';
        if i > 1
            Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';
            objV = objV + sum(shat(:));
        end
    end 
 end

Y = ifft(Yhat,[],3);  %% 逆傅里叶变换
if mode == 1
    X = Yi2X(Y,3);
elseif mode == 3
    X = shiftdim(Y, 2);  %% 调整张量维度[n,n,k]
else
    X = Y;
end

x = X(:);  %% 输出为列向量

end
 