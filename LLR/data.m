clear;
% rng(1000);
L = 500; M = 20; %screen 1
mag = 3; nGWAS = 20; OddsRatio = 1/49;
n = 5000; LDfolder = 'tmp1';
rho = 0.5; h2 = 0.5;

wd = 'D:\Documents\LLR' 
% 'G:\Data\GoogleDrive\Work\GPA_matlab\OtherWork\LowRank\Z_score\Duke-NUS\Allfiles\packages_llr';%wd = '/home/svu/gmsjl/LLR/LLR2';
nstep = 3000; eps = 0.1;nref = 400;
% wdPAINTOR = '~/Work/LLR';
cd(wd);
maf = load('maf_10000.txt');
nsnp = L*M;
rng(1000);

[TrueSnp,Zstat,loci,y,~,RR] = generateData3(maf,L,M,rho,h2,n,nref, mag,OddsRatio,nGWAS);