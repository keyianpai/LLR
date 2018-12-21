function [TrueSnp,Zstat,loci,y,SIGMA,RR] = generateData3(maf,L,M,rho,h2,n,nref, mag,OddsRatio,nGWAS)

SIGMA = zeros(M,M);
for j = 1:M
    for k= 1:M
        SIGMA(j,k) = rho^(abs(j-k));
    end
end

r = 2;
c1=[1;-1];%（2,1）列向量,分号，换行
c2=[1;1];
v = [c1*ones(1,nGWAS/2), c2*ones(1,nGWAS/2)];%GWAS研究的表型分成两组
x = randn(L,r)*v;%x：low rank matrix
x0 = ones(L,nGWAS)*log(OddsRatio);
prob = 1./(1+exp(-x*mag-x0));%矩阵中每个元素可以单独来理解，prob矩阵表示L（500）个位点与第k种表型有联系的概率

y = binornd(1,prob);%的第k列表示L（500）个位点与第k种表型是否有联系，有联系为1，按prob中的概率产生

loci = ones(M,1)*(1:L);
loci=loci(:);%20个1,20个2，。。。，20个500

nsnp = L*M;
TrueSnp = zeros(nsnp, nGWAS);

for k = 1:nGWAS
    for l = 1: L
        if y(l,k) == 1
            locusIndx = zeros(M,1);
            index = (M*(l-1)+1): (M*l);
            indx = randsample(M,1);
            locusIndx(indx) = 1;
            TrueSnp(index,k) = locusIndx;
        end
    end
end
%每次循环，若第l个位点与第k个表型有联系，即y(l,k) ==1，则将此位点的M个snp中随机抽一个作为与此表型有联系的snp，循环结束生成TrueSnp矩阵

b = zeros(nsnp,nGWAS);
for k = 1:nGWAS
    b(TrueSnp(:,k)==1,k) = randn(sum(TrueSnp(:,k)==1),1);
end
%b矩阵与TrueSnp矩阵一样大，TrueSnp为1的位置生成一个标准正态分布的随机数，
h2_true = zeros(nGWAS,1);
betah = zeros(nsnp,nGWAS); %Zstat = betah./sqrt(s2);
s2 = zeros(nsnp,nGWAS);
% pvals = zeros(nsnp,nGWAS);

X00 = cell(L,1); Y  = cell(nGWAS,1); % Xo = cell(nGWAS,1);
% RRaw = cell(L,2); 
stderror = zeros(nGWAS,1);

for k = 1:nGWAS
    Y0 = 0;
    %Xo{k} = [];
    for l = 1:L
        index = (M*(l-1)+1): (M*l);
        AAprob = maf(index).^2.;%AA count as 2,AAprob(1,M)
        Aaprob = 2*maf(index).*(1-maf(index));
        quanti = [1-Aaprob-AAprob; 1- AAprob]';%(M,2)
        X00{l} = mvnrnd(zeros(M,1),SIGMA,n);%对每个性状而言这个循环的内容是一样的，似乎与k无关，唯一不同在于产生5000个M维（一个位点内snp个数）正态分布向量，作为5000个个体在l位点处基因型的模拟，位点内snp关系体现在SIGMA矩阵
        Xrefc = zeros(n,M);
        for j = 1:M
            cutoff = norminv(quanti(j,:));%(1,2)，Inverse of the normal cumulative distribution function (cdf).
            Xrefc(X00{l}(:,j) < cutoff(1),j) = 0;
            Xrefc(X00{l}(:,j) >= cutoff(1) & X00{l}(:,j) < cutoff(2),j) = 1;
            Xrefc(X00{l}(:,j) >= cutoff(2),j) = 2;
        end
        %内层循环对X00矩阵进行离散化，使得AA，Aa，aa频率满足一定要求，A为maf对应基因型，离散化后的0,1,2,指的是A的数量
        X00{l} = Xrefc;
        X = X00{l};%(X00{l} - repmat(mean(X00{l}),n,1))./repmat(std(X00{l}),n,1)/sqrt((alpha)*nsnp);
        Y0 = Y0 + X*b(loci==l,k);%X：（5000，M）
        % RRaw{l,k} = corr(X00{l});
        % Xo{k} = [Xo{k},Xrefc];
    end
    stderror(k) = std(Y0)*sqrt((1-h2)/h2);
    e = stderror(k)*randn(n,1);%/mag(k);
    Y{k} = Y0 + e;
    %以b矩阵为系数，与离散化的X矩阵（5000，M)相乘得到5000个个体的性状Y0，Y0加上噪声得到Y{k}
    %上下两个1：L循环，上面为数据生成过程，下面为参数估计过程
    h2_true(k) = var(Y0)/var(Y{k});
    for l = 1:L
        tmpbeta = zeros(M,1); tmps2 = zeros(M,1); tmppvals = zeros(M,1);
        for i = 1:M
            [betahat,se2,~,pval,~] = LinearRegCan2(X00{l}(:,i),Y{k});
            tmpbeta(i) = betahat(1);
            tmps2(i) = se2(1);
            tmppvals(i) = pval(1);
        end
        betah(loci==l,k) = tmpbeta;
        s2(loci==l,k) = tmps2;
        % pvals(loci==l,k) = tmppvals;
    end    
end

clear X00 Y stderror Xrefc X;

R = []; RR = cell(L,1);
for l = 1:L
    index = (M*(l-1)+1): (M*l);
    AAprob = maf(index).^2.;
    Aaprob = 2*maf(index).*(1-maf(index));
    quanti = [1-Aaprob-AAprob; 1- AAprob]';
    Xref = mvnrnd(zeros(M,1),SIGMA,nref);
    Xrefc = zeros(n,M);
    for j = 1:M
        cutoff = norminv(quanti(j,:));
        Xrefc(Xref(:,j) < cutoff(1),j) = 0;
        Xrefc(Xref(:,j) >= cutoff(1) & Xref(:,j) < cutoff(2),j) = 1;
        Xrefc(Xref(:,j) >= cutoff(2),j) = 2;
    end
    RR{l} = corr(Xrefc);
    % R = blkdiag(R,RR{l}); 
end

Zstat = betah./sqrt(s2);
clear Xref Xrefc;