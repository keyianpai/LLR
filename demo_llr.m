clear;
% rng(1000);
L = 500; M = 20; %screen 1 L：位点个数，M：每个位点内snp的个数
mag = 3; nGWAS = 20; OddsRatio = 1/49;%mag：由X矩阵和x0向量产生prob矩阵时X前的系数，（由prob继而产生y），nGWAS;表型研究的个数，log(OddsRatio)为x0的系数
n = 5000; LDfolder = 'tmp1';%n：模拟人群个数
rho = 0.5; h2 = 0.5;%rho决定位点snp的协方差矩阵SIGMA

wd = 'D:\Documents\LLR'; 
% 'G:\Data\GoogleDrive\Work\GPA_matlab\OtherWork\LowRank\Z_score\Duke-NUS\Allfiles\packages_llr';%wd = '/home/svu/gmsjl/LLR/LLR2';
nstep = 3000; eps = 0.1;nref = 400;%nstep：EM算法迭代步数，eps：学习率，nref：为了估计每个位点snp的协方差矩阵，模拟产生的参考人群个数
% wdPAINTOR = '~/Work/LLR';
cd(wd);
maf = load('maf_10000.txt');%10000个snp的maf，论文中说由均匀分布产生
nsnp = L*M;%nsnp：snp的总个数，L：位点个数，M：每个位点内snp的个数
rng(1000);

[TrueSnp,Zstat,loci,y,~,RR] = generateData3(maf,L,M,rho,h2,n,nref, mag,OddsRatio,nGWAS);% y is produced by X
%TrueSnp：表示snp位点是否与表型有联系的0,1矩阵，1表示有联系，有nsnp行，nGWAS列
Cpio = sum(y,1)/L;%y 的第k列表示L（500）个位点与第k种表型是否有联系（由prob矩阵转化而来，prob矩阵表示L（500）个位点与第k种表型有联系的概率），每一列求和再除以每一列的元素个数L（也就是行数）得到nGWAS（20）个表型与所有位点有联系的平均概率。
    
LD = cell(L,1);
for j = 1:L
    LD{j} = RR{j};
end

%%LLR
options = GPA_lowrank_set(nGWAS,[]);
options.eps = eps;
options.maxIters = nstep;
options.ncp = 1;

obj0 = Zscore_init(Zstat, loci, LD, 3.7,options);
Cpi1 = obj0.Cpi1;

nfold =5;
[loglik_masked, ave_crit] = cv_Zscore_lowrank_boosting(Zstat, LD, 3.7, loci, obj0, options, nfold);
save('ave_crit.txt','ave_crit','-ascii','-tabs');
save('loglik_masked.txt','loglik_masked','-ascii','-tabs');

[~,ind] = max(ave_crit);
index = 1:nstep;
optstep = min(index(ave_crit > ave_crit(ind) - 5));

options.maxIters = optstep; mask = ones(L,nGWAS);
obj = Zscore_lowrank(Zstat, loci, LD, 3.7, obj0, mask,options);
x = obj.x;
[~,~,V] = svd(x, 'econ');

subplot(2,2,1)
plot(obj.Dall')
plot(ave_crit)
scatter(V(1,:),V(:,2));
Traits = num2cell(1:20);
text(V(1,:)+0.02, V(:,2)+0.02, Traits);

d=pdist(V(:,1:2));
z=linkage(d);
dendrogram(z);
