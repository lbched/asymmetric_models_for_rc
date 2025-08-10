%% In-Sample Estimation Results & Hypothesis Testing (Tables:3-9;16-24)
clear
clc
warning('off','all');

%% The file paths are relative to the folder Data/
% Loading Data
% Realized Covariance series
rc = readmatrix('RC.csv')*25200; % annualized in percentage

T = length(rc); % sample size
n = 6;  % cross-sectional dimension

RC = zeros(n,n,T);
for t = 1:T
   RC(:,:,t) = buildSymmetric(rc(t,:),n);
end

% Close-to-Close return-based components (via eq.3)
pc = readmatrix('PC.csv')*25200;
nc = readmatrix('NC.csv')*25200;
mc = readmatrix('MC.csv')*25200;
mpc = readmatrix('MPC.csv')*25200;
mnc = readmatrix('MNC.csv')*25200; 

pmc = zeros(T,n*(n+1)/2);
for t = 1:T
   pmc(t,:) = pc(t,:)+mc(t,:);
end

PMC = zeros(n,n,T);
for t = 1:T
   PMC(:,:,t) = buildSymmetric(pmc(t,:),n);
end

PC = zeros(n,n,T);
for t = 1:T
   PC(:,:,t) = buildSymmetric(pc(t,:),n);
end

NC = zeros(n,n,T);
for t = 1:T
   NC(:,:,t) = buildSymmetric(nc(t,:),n);
end

MC = zeros(n,n,T);
for t = 1:T
   MC(:,:,t) = buildSymmetric(mc(t,:),n);
end

MPC = zeros(n,n,T);
for t = 1:T
   MPC(:,:,t) = buildSymmetric(mpc(t,:),n);
end

MNC = zeros(n,n,T);
for t = 1:T
   MNC(:,:,t) = buildSymmetric(mnc(t,:),n);
end

% Open-to-Close return-based components (via eq.3)
pc_oc = readmatrix('PC_oc.csv')*25200;
nc_oc = readmatrix('NC_oc.csv')*25200;
mc_oc = readmatrix('MC_oc.csv')*25200;
mpc_oc = readmatrix('MPC_oc.csv')*25200;
mnc_oc = readmatrix('MNC_oc.csv')*25200;

pmc_oc = zeros(T,n*(n+1)/2);
for t = 1:T
   pmc_oc(t,:) = pc_oc(t,:)+mc_oc(t,:);
end

PMC_oc = zeros(n,n,T);
for t = 1:T
   PMC_oc(:,:,t) = buildSymmetric(pmc_oc(t,:),n);
end

PC_oc = zeros(n,n,T);
for t = 1:T
   PC_oc(:,:,t) = buildSymmetric(pc_oc(t,:),n);
end 

NC_oc = zeros(n,n,T);
for t = 1:T
   NC_oc(:,:,t) = buildSymmetric(nc_oc(t,:),n);
end 

MC_oc = zeros(n,n,T);
for t = 1:T
   MC_oc(:,:,t) = buildSymmetric(mc_oc(t,:),n);
end

MPC_oc = zeros(n,n,T);
for t = 1:T
   MPC_oc(:,:,t) = buildSymmetric(mpc_oc(t,:),n);
end

MNC_oc = zeros(n,n,T);
for t = 1:T
   MNC_oc(:,:,t) = buildSymmetric(mnc_oc(t,:),n);
end

% Semi-covariance components (via eq. 10)
psc = readmatrix('PSC.csv')*25200;
nsc = readmatrix('NSC.csv')*25200;
msc = readmatrix('MSC.csv')*25200;
mpsc = readmatrix('MPSC.csv')*25200;
mnsc = readmatrix('MNSC.csv')*25200;

PSC = zeros(n,n,T);
for t = 1:T
   PSC(:,:,t) = buildSymmetric(psc(t,:),n);
end

NSC = zeros(n,n,T);
for t = 1:T
   NSC(:,:,t) = buildSymmetric(nsc(t,:),n);
end

MSC = zeros(n,n,T);
for t = 1:T
   MSC(:,:,t) = buildSymmetric(msc(t,:),n);
end

MPSC = zeros(n,n,T);
for t = 1:T
   MPSC(:,:,t) = buildSymmetric(mpsc(t,:),n);
end

MNSC = zeros(n,n,T);
for t = 1:T
   MNSC(:,:,t) = buildSymmetric(mnsc(t,:),n);
end

%% The file paths are relative to the folder Fun/
%% Main output of the estimations, function "Max_lik": 
% Coefficient values: beta_(d/p)_name_of_the_model
% LLF value: logl_(d/p)_name_of_the_model
% AIC value: aic_(d/p)_name_of_the_model
% BIC value: bic_(d/p)_name_of_the_model
% % E.g., fitted coefficients betasym for scalar, betadsym for diagonal,
% % betapsym for PLT version of the sym model

%% Scalar models 
% Table 16 (CC-return and semi-decomposition-based specifications);
% Table 19 (OC-return-based specifications); 
% Table 24 (extended asymmetric specification); 
% Table 3 (rows: LLF-scalar, AIC-scalar, BIC-scalar)
% sym
lb = zeros(2,1)';
ub = ones(2,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Scsymllk(p,RC(:,:,1:end),n,T);
x0 = [0.4; 0.6]';
[betasym,stderrsym,vcsym,loglsym] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aicsym,bicsym] = aicbic(loglsym,2,T,Normalize=true);
% Table 16 (column sym): coefficients: betasym 
% Table 3&16 (column sym): LLF-scalar: loglsym
% Table 3&16 (column sym): AIC-scalar: aicsym
% Table 3&16 (column sym): BIC-scalar: bicsym

% tr 
lb = zeros(3,1)';
ub = ones(3,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Sctrllk(p,RC(:,:,1:end),PMC(:,:,1:end),NC(:,:,1:end),n,T);
x0 = [0.1; 0.07; 0.8]';
[betatr,stderrtr,vctr,logltr] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aictr,bictr] = aicbic(logltr,3,T,Normalize=true);
% Table 16 (column tr): coefficients: betatr 
% Table 3&16 (column tr): LLF-scalar: logltr
% Table 3&16 (column tr): AIC-scalar: aictr
% Table 3&16 (column tr): BIC-scalar: bictr

% trPNM
lb=[0;0;-1;0]';
ub=ones(4,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) SctrPNMllk(p,RC(:,:,1:end),PC(:,:,1:end),NC(:,:,1:end),MC(:,:,1:end),n,T);
x0 = [0.11;0.17;0.13;0.8]';
[betatrPNM,stderrtrPNM,vctrPNM,logltrPNM]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aictrPNM,bictrPNM] = aicbic(logltrPNM,4,T,Normalize=true);

% trPNtauM
lb=[0;0;-1;-1;0]';
ub=ones(5,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) SctrPNtauMllk(p,RC(:,:,1:end),PC(:,:,1:end),NC(:,:,1:end),MPC(:,:,1:end),MNC(:,:,1:end),n,T);
x0 = [0.11;0.19;0.13;0.1;0.8]';
[betatrPNtauM,stderrtrPNtauM,vctrPNtauM,logltrPNtauM]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aictrPNtauM,bictrPNtauM] = aicbic(logltrPNtauM,5,T,Normalize=true);

% semi
lb=[0;0;-1;0]';
ub=ones(4,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Scsemillk(p,RC(:,:,1:end),PSC(:,:,1:end),NSC(:,:,1:end),MSC(:,:,1:end),n,T);
x0 = [0.21;0.27;0.1;0.8]';
[betasemi,stderrsemi,vcsemi,loglsemi]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aicsemi,bicsemi] = aicbic(loglsemi,4,T,Normalize=true);

% semi-tau
lb=[0;0;-1;-1;0]';
ub=ones(5,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Scsemitaullk(p,RC(:,:,1:end),PSC(:,:,1:end),NSC(:,:,1:end),MPSC(:,:,1:end),MNSC(:,:,1:end),n,T);
x0 = [0.21;0.27;0.05;0.03;0.8]';
[betasemitau,stderrsemitau,vcsemitau,loglsemitau]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aicsemitau,bicsemitau] = aicbic(loglsemitau,5,T,Normalize=true);

% tr^oc
lb = zeros(3,1)';
ub = ones(3,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Sctrllk(p,RC(:,:,1:end),PMC_oc(:,:,1:end),NC_oc(:,:,1:end),n,T);
x0 = [0.1; 0.07; 0.8]';
[betatr_oc,stderrtr_oc,vctr_oc,logltr_oc] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aictr_oc,bictr_oc] = aicbic(logltr_oc,3,T,Normalize=true);
% Table 19 (column tr): coefficients: betatr_oc
% Table 3&19 (column tr(^oc)):: LLF-scalar: logltr_oc
% Table 3&19 (column tr(^oc)): AIC-scalar: aictr_oc
% Table 3&19 (column tr(^oc)): BIC-scalar: bictr_oc

% trPNM^oc
lb=[0;0;-1;0]';
ub=ones(4,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) SctrPNMllk(p,RC(:,:,1:end),PC_oc(:,:,1:end),NC_oc(:,:,1:end),MC_oc(:,:,1:end),n,T);
x0 = [0.11;0.17;0.13;0.8]';
[betatrPNM_oc,stderrtrPNM_oc,vctrPNM_oc,logltrPNM_oc]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aictrPNM_oc,bictrPNM_oc] = aicbic(logltrPNM_oc,4,T,Normalize=true);
% Table 19 (column trPNM): coefficients: betatrPNM_oc
% Table 3&19 (column trPNM(^oc)): LLF-scalar: logltrPNM_oc
% Table 3&19 (column trPNM(^oc)): AIC-scalar: aictrPNM_oc
% Table 3&19 (column trPNM(^oc)): BIC-scalar: bictrPNM_oc

% trPNtauM^oc
lb=[0;0;-1;-1;0]';
ub=ones(5,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) SctrPNtauMllk(p,RC(:,:,1:end),PC_oc(:,:,1:end),NC_oc(:,:,1:end),MPC_oc(:,:,1:end),MNC_oc(:,:,1:end),n,T);
x0 = [0.11;0.19;0.13;0.1;0.8]';
[betatrPNtauM_oc,stderrtrPNtauM_oc,vctrPNtauM_oc,logltrPNtauM_oc]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aictrPNtauM_oc,bictrPNtauM_oc] = aicbic(logltrPNtauM_oc,5,T,Normalize=true);

% trPNM-semi
lb=[0;0;-1;0;0;-1;0]';
ub=ones(7,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) SctrPNM_semillk(p,RC(:,:,1:end),PC(:,:,1:end),NC(:,:,1:end),MC(:,:,1:end),PSC(:,:,1:end),NSC(:,:,1:end),MSC(:,:,1:end),n,T);
x0 = [0.11;0.17;0.13;0.21;0.27;0.1;0.8]';
[betatrPNM_semi,stderrtrPNM_semi,vctrPNM_semi,logltrPNM_semi]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aictrPNM_semi,bictrPNM_semi] = aicbic(logltrPNM_semi,7,2517,Normalize=true);
% Table 24 (column trPNM-semiM): coefficients: betatrPNM_semi
% Table 24 (column trPNM-semi): LLF: logltrPNM_semi
% Table 24 (column trPNM-semi): AIC: aictrPNM_semi
% Table 24 (column trPNM-semi): BIC: bictrPNM_semi

% trPNtauM-semi-tau
lb=[0;0;-1;-1;0;0;-1;-1;0]';
ub=ones(9,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) SctrPNtauM_semitaullk(p,RC(:,:,1:end),PC(:,:,1:end),NC(:,:,1:end),MPC(:,:,1:end),MNC(:,:,1:end),PSC(:,:,1:end),NSC(:,:,1:end),MPSC(:,:,1:end),MNSC(:,:,1:end),n,T);
x0 = [0.11;0.19;0.13;0.1;0.21;0.27;0.05;0.03;0.8]';
[betatrPNtauM_semitau,stderrtrPNtauM_semitau,vctrPNtauM_semitau,logltrPNtauM_semitau]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aictrPNtauM_semitau,bictrPNtauM_semitau] = aicbic(logltrPNtauM_semitau,9,2517,Normalize=true);
% Table 24 (column trPNtauM-semi-tau): coefficients: betatrPNtauM_semitau
% Table 24 (column trPNtauM-semi-tau): LLF: logltrPNtauM_semitau
% Table 24 (column trPNtauM-semi-tau): AIC: aictrPNtauM_semitau
% Table 24 (column trPNtauM-semi-tau): BIC: bictrPNtauM_semitau

%% Selected Diagonal models
% Table 17 (CC-return and semi-decomposition-based specifications);
% Table 20 (OC-return-based specifications); 
% Table 3 (rows: LLF-diagonal, AIC-diagonal, BIC-diagonal)
% sym
lb=zeros(n*2,1)';
ub=ones(n*2,1)';
options=optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk= @(p) Diagsymllk(p,RC(:,:,1:end),n,T);
x0=[sqrt(0.27)*ones(n,1); sqrt(0.7)*ones(n,1)]';
[betadsym,stderrdsym,vcdsym,logldsym]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aicdsym,bicdsym] = aicbic(logldsym,12,T,Normalize=true);
stsdsym = zeros(12,1);
for t = 1:12
   stsdsym(t) = betadsym(t)/(stderrdsym(t));
end
% Table 17 (column sym): coefficients: betadsym 
% Table 17 (column sym): coefficient significance: stsdsym
% Table 3&17 (column sym): LLF-diagonal: logldsym
% Table 3&17 (column sym): AIC-diagonal: aicdsym
% Table 3&17 (column sym): BIC-diagonal: bicdsym

% tr
lb=zeros(n*3,1)';
ub=ones(n*3,1)';
options=optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk= @(p) Diagtrllk(p,RC(:,:,1:end),PMC(:,:,1:end),NC(:,:,1:end),n,T);
x0=[sqrt(0.24)*ones(n,1); sqrt(0.28)*ones(n,1); sqrt(0.71)*ones(n,1)]';
[betadtr,stderrdtr,vcdtr,logldtr]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aicdtr,bicdtr] = aicbic(logldtr,18,T,Normalize=true);
stsdtr = zeros(18,1);
for t = 1:18
   stsdtr(t) = betadtr(t)/(stderrdtr(t));
end
% Table 17 (column tr): coefficients: betadtr
% Table 17 (column tr): coefficient significance: stsdtr
% Table 3&17 (column tr): LLF-diagonal: logldtr
% Table 3&17 (column tr): AIC-diagonal: aicdtr
% Table 3&17 (column tr): BIC-diagonal: bicdtr

% tr^oc
lb=zeros(n*3,1)';
ub=ones(n*3,1)';
options=optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk= @(p) Diagtrllk(p,RC(:,:,1:end),PMC_oc(:,:,1:end),NC_oc(:,:,1:end),n,T);
x0=[sqrt(0.24)*ones(n,1); sqrt(0.28)*ones(n,1); sqrt(0.71)*ones(n,1)]';
[betadtr_oc,stderrdtr_oc,vcdtr_oc,logldtr_oc]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aicdtr_oc,bicdtr_oc] = aicbic(logldtr_oc,18,T,Normalize=true);
stsdtr_oc = zeros(18,1);
for t = 1:18
   stsdtr_oc(t) = betadtr_oc(t)/(stderrdtr_oc(t));
end
% Table 20 (column tr): coefficients: betadtr_oc
% Table 20 (column tr): coefficient significance: stsdtr_oc
% Table 3&20 (column tr(^oc)): LLF-diagonal: logldtr_oc
% Table 3&20 (column tr(^oc)): AIC-diagonal: aicdtr_oc
% Table 3&20 (column tr(^oc)): BIC-diagonal: bicdtr_oc

%% Selected PLT models
% Table 18 (CC-return and semi-decomposition-based specifications);
% Table 21 (OC-return-based specifications); 
% Table 3 (rows: LLF-PLT, AIC-PLT, BIC-PLT)
% sym
lb=[0 -0.3 -0.3 -0.3 -0.3 -0.3 0 0 0 0 0 0 0 0 0 0 0]';
ub=ones(17,1)';
options=optimset('Display','iter','MaxFunEvals',1000000,'MaxIter',100000);
llk= @(p) PLTsymllk(p,RC(:,:,1:end),n,T);
x0=[0.13;0.02;0.02;0.02;0.02;0.02;0.15;0.15;0.15;0.15;0.15;0.5;0.5;0.5;0.5;0.5;0.5]';
[betapsym,stderrpsym,vcpsym,loglpsym]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aicpsym,bicpsym] = aicbic(loglpsym,17,T,Normalize=true);
stspsym = zeros(17,1);
for t = 1:17
   stspsym(t) = betapsym(t)/(stderrpsym(t));
end
% Table 18 (column sym): coefficients: betapsym 
% Table 18 (column sym): coefficient significance: stspsym
% Table 3&18 (column sym): LLF-PLT: loglpsym
% Table 3&18 (column sym): AIC-PLT: aicpsym
% Table 3&18 (column sym): BIC-PLT: bicpsym

% tr 
lb=[0 -0.05 -0.05 -0.05 -0.05 -0.05 0 0 0 0 0 0 -0.1 -0.1 -0.1 -0.1 -0.1 0 0 0 0 0 0 0 0 0 0 0]';
ub=ones(28,1)';
options=optimset('Display','iter','MaxFunEvals',1000000,'MaxIter',100000);
llk= @(p) PLTtrllk(p,RC(:,:,1:end),PMC(:,:,1:end),NC(:,:,1:end),n,T);
x0=[0.13;0.02;0.02;0.02;0.02;0.02;0.15;0.15;0.15;0.15;0.15;0.19;0.01;0.01;0.01;0.01;0.01;0.19;0.19;0.19;0.19;0.19;0.5;0.5;0.5;0.5;0.5;0.5]';
[betaptr,stderrptr,vcptr,loglptr]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aicptr,bicptr] = aicbic(loglptr,28,T,Normalize=true);
stsptr = zeros(28,1);
for t = 1:28
   stsptr(t) = betaptr(t)/(stderrptr(t));
end
% Table 18 (column tr): coefficients: betaptr
% Table 18 (column tr): coefficient significance: stsptr
% Table 3&18 (column tr): LLF-PLT: loglptr
% Table 3&18 (column tr): AIC-PLT: aicptr
% Table 3&18 (column tr): BIC-PLT: bicptr

% tr^oc
lb=[0 -0.1 -0.1 -0.1 -0.1 -0.1 0 0 0 0 0 0 -0.1 -0.1 -0.1 -0.1 -0.1 0 0 0 0 0 0 0 0 0 0 0]';
ub=ones(28,1)';
options=optimset('Display','iter','MaxFunEvals',1000000,'MaxIter',100000);
llk= @(p) PLTtrllk(p,RC(:,:,1:end),PMC_oc(:,:,1:end),NC_oc(:,:,1:end),n,T);
x0=[0.368804558362593;0.01;0.01;0.01;0.01;0.01;0.538065952083568;0.530318973022856;0.511549377512276;0.545485168848220;0.573257572781041;0.472004744259667;0.01;0.01;0.01;0.01;0.01;0.567506080656019;0.554722139862587;0.533979918990755;0.577588629369931;0.611362715636051;0.896864395332714;0.791443957839432;0.809195896886628;0.815725618142746;0.786526806403777;0.766984572255703]';
[betaptr_oc,stderrptr_oc,vcptr_oc,loglptr_oc]=Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
[aicptr_oc,bicptr_oc] = aicbic(loglptr_oc,28,T,Normalize=true);
stsptr_oc = zeros(28,1);
for t = 1:28
   stsptr_oc(t) = betaptr_oc(t)/(stderrptr_oc(t));
end
% Table 21 (column tr): coefficients: betaptr_oc
% Table 21 (column tr): coefficient significance: stsptr_oc
% Table 3&21 (column tr(^oc)): LLF-PLT: loglptr_oc
% Table 3&21 (column tr(^oc)): AIC-PLT: aicptr_oc
% Table 3&21 (column tr(^oc)): BIC-PLT: bicptr_oc

%% Table 4: Selected LR Tests
% Scalar sym vs. tr (row 3, column 4)
[h_sc_sym_tr,p_sc_sym_tr,stat_sc_sym_tr]=lratiotest(logltr,loglsym,1);
% Scalar sym vs. tr^oc (row 4, column 4)
[h_sc_sym_tr_oc,p_sc_sym_tr_oc,stat_sc_sym_tr_oc]=lratiotest(logltr_oc,loglsym,1);
% Diagonal sym vs. tr (row 3, column 5)
[h_d_sym_tr,p_d_sym_tr,stat_d_sym_tr]=lratiotest(logldtr,logldsym,6);
% Diagonal sym vs. tr^oc (row 4, column 5)
[h_d_sym_tr_oc,p_d_sym_tr_oc,stat_d_sym_tr_oc]=lratiotest(logldtr_oc,logldsym,6);
% PLT sym vs. tr (row 3, column 6)
[h_plt_sym_tr,p_plt_sym_tr,stat_plt_sym_tr]=lratiotest(loglptr,loglpsym,11);
% PLT sym vs. tr^oc (row 4, column 6)
[h_plt_sym_tr_oc,p_plt_sym_tr_oc,stat_plt_sym_tr_oc]=lratiotest(loglptr_oc,loglpsym,11);

%% Table 5: Selected LR Tests
% Scalar vs. diagonal sym (column sym, rows 2-4)
[h_sc_d_sym,p_sc_d_sym,stat_sc_d_sym]=lratiotest(logldsym,loglsym,10);
% Scalar vs. PLT sym (column sym, rows 5-7)
[h_sc_plt_sym,p_sc_plt_sym,stat_sc_plt_sym]=lratiotest(loglpsym,loglsym,15);
% Diagonal vs. PLT sym (column sym, rows 8-10)
[h_d_plt_sym,p_d_plt_sym,stat_d_plt_sym]=lratiotest(loglpsym,logldsym,5);
% Scalar vs. diagonal tr (column tr, rows 2-4)
[h_sc_d_tr,p_sc_d_tr,stat_sc_d_tr]=lratiotest(logldtr,logltr,15);
% Scalar vs. PLT tr (column tr, rows 5-7)
[h_sc_plt_tr,p_sc_plt_tr,stat_sc_plt_tr]=lratiotest(loglptr,logltr,25);
% Diagonal vs. PLT tr (column tr, rows 8-10)
[h_d_plt_tr,p_d_plt_tr,stat_d_plt_tr]=lratiotest(loglptr,logldtr,10);
% Scalar vs. diagonal tr_oc (column tr^oc, rows 2-4)
[h_sc_d_tr_oc,p_sc_d_tr_oc,stat_sc_d_tr_oc]=lratiotest(logldtr_oc,logltr_oc,15);
% Scalar vs. PLT tr_oc (column tr^oc, rows 5-7)
[h_sc_plt_tr_oc,p_sc_plt_tr_oc,stat_sc_plt_tr_oc]=lratiotest(loglptr_oc,logltr_oc,25);
% Diagonal vs. PLT tr_oc (column tr^oc, rows 8-10)
[h_d_plt_tr_oc,p_d_plt_tr_oc,stat_d_plt_tr_oc]=lratiotest(loglptr_oc,logldtr_oc,10);

%% Table 6: Variance Coefficients, p.1
% Scalar sym (row 3)
ssym_coeff=betasym.^2;
% Diagonal sym 
dsym_coeff=betadsym.^2;
% Asset 1 (row 4)
dsym_coeff(1)
dsym_coeff(7)
% Assets 2-6 (row 5)
mean(dsym_coeff(2:6))
mean(dsym_coeff(8:12))
% PLT sym
psym_coeff=betapsym.^2;
% Asset 1 (row 6)
psym_coeff(1)
psym_coeff(12)
% Assets 2-6 (row 7)
mean(psym_coeff(7:11))
mean(psym_coeff(13:17))
% Scalar tr (row 8)
str_coeff=betatr.^2;
% Diagonal tr
dtr_coeff=betadtr.^2;
% Asset 1 (row 9)
dtr_coeff(1)
dtr_coeff(7)
dtr_coeff(13)
% Assets 2-6 (row 10)
mean(dtr_coeff(2:6))
mean(dtr_coeff(8:12))
mean(dtr_coeff(14:18))
% PLT tr
ptr_coeff=betaptr.^2;
% Asset 1 (row 11)
ptr_coeff(1)
ptr_coeff(12)
ptr_coeff(23)
% Assets 2-6 (row 12)
mean(ptr_coeff(7:11))
mean(ptr_coeff(18:22))
mean(ptr_coeff(24:28))

%% Table 7: Variance Coefficients, p.2
% PLT sym 
psym_coeff=betapsym.^2;
% Asset 2-6; positive & negative component (column sym, rows 2-3)
mean(psym_coeff(2:6))
% Asset 2-6; positive & negative component (column sym, rows 4-6)
mean(2*betapsym(2:6).*betapsym(7:11))
% PLT tr
ptr_coeff=betaptr.^2;
% Asset 2-6; positive component (column tr, row 2)
mean(ptr_coeff(2:6)) 
% Assets 2-6; negative component (column tr, row 3)
mean(ptr_coeff(13:17))
% Asset 2-6; positive component (column tr, row 4)
mean(2*betaptr(2:6).*betaptr(7:11))
% Asset 2-6; negative component (column tr, row 5)
mean(2*betaptr(13:17).*betaptr(18:22))

%% Table 8: Covariance Coefficients, p.1
% Scalar sym (row 3)
ssym_coeff_cov(1)=betasym(1)^2;
% Diagonal sym
dsym_coeff_cov=betadsym;
% Assets 1-other (row 4)
mean(dsym_coeff_cov(1).*dsym_coeff_cov(2:6))
% Assets 2-6 (row 5)
mean(dsym_coeff_cov(2).*dsym_coeff_cov(3:6))
% PLT sym
psym_coeff_cov=betapsym;
% Assets 1-other (row 6)
mean(psym_coeff_cov(1).*psym_coeff_cov(7:11))
% Assets 2-6 (row 7)
mean(psym_coeff_cov(7).*psym_coeff_cov(8:11))
% Scalar tr (row 8)
str_coeff_cov(1:2)=betatr(1:2).^2;
% Diagonal tr
dtr_coeff_cov=betadtr;
% Assets 1-other (row 9)
mean(dtr_coeff_cov(1).*dtr_coeff_cov(2:6))
mean(dtr_coeff_cov(7).*dtr_coeff_cov(8:12))
% Assets 2-6 (row 10)
mean(dtr_coeff_cov(2).*dtr_coeff_cov(3:6))
mean(dtr_coeff_cov(8).*dtr_coeff_cov(9:12))
% PLT tr
ptr_coeff_cov=betaptr;
% Assets 1-other (row 11)
mean(ptr_coeff_cov(1).*ptr_coeff_cov(7:11))
mean(ptr_coeff_cov(12).*ptr_coeff_cov(18:22))
% Assets 2-6 (row 12)
mean(ptr_coeff_cov(7).*ptr_coeff_cov(8:11))
mean(ptr_coeff_cov(18).*ptr_coeff_cov(19:22))

%% Table 9: Covariance Coefficients, p.2
% Scalar sym (row 3)
ssym_coeff_cov_b=betasym(2)^2;
% Diagonal sym
dsym_coeff_cov=betadsym;
% Assets 1-other (row 4)
mean(dsym_coeff_cov(7).*dsym_coeff_cov(8:12))
% Assets 2-6 (row 5)
mean(dsym_coeff_cov(8).*dsym_coeff_cov(9:12))
% PLT sym
psym_coeff_cov=betapsym;
% Assets 1-other (row 6)
mean(psym_coeff_cov(12).*psym_coeff_cov(13:17))
mean(psym_coeff_cov(1).*psym_coeff_cov(2:6))
% Assets 2-6 (row 7)
mean(psym_coeff_cov(13).*psym_coeff_cov(14:17))
% Scalar tr (row 8)
str_coeff_cov_beta=betatr(3)^2;
% Diagonal tr
dtr_coeff_cov=betadtr;
% Assets 1-other (row 9)
mean(dtr_coeff_cov(13).*dtr_coeff_cov(14:18))
% Assets 2-6 (row 10)
mean(dtr_coeff_cov(14).*dtr_coeff_cov(15:18))
% PLT tr
ptr_coeff_cov=betaptr;
% Assets 1-other (row 11)
mean(ptr_coeff_cov(23).*ptr_coeff_cov(24:28))
mean(ptr_coeff_cov(1).*ptr_coeff_cov(2:6))
mean(ptr_coeff_cov(12).*ptr_coeff_cov(13:17))
% Assets 2-6 (row 12)
mean(ptr_coeff_cov(24).*ptr_coeff_cov(25:28))

%% Table 22: Leverage effects in variance equations
% Scalar tr
% Test statistics Variance coefficients (column 5) 
g=[-2*betatr(1),2*betatr(2),0]';
dcov=g'*vctr*g;
ts_s=(betatr(2)^2-betatr(1)^2)/sqrt(dcov);
ts_s;
% Diagonal tr
% Test statistics Variance coefficients (column 5) 
g=[-2*betadtr(1),0,0,0,0,0,2*betadtr(7),0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d1=(betadtr(7)^2-betadtr(1)^2)/sqrt(dcov);
ts_d1;

g=[0,-2*betadtr(2),0,0,0,0,0,2*betadtr(8),0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d2=(betadtr(8)^2-betadtr(2)^2)/sqrt(dcov);
ts_d2;

g=[0,0,-2*betadtr(3),0,0,0,0,0,2*betadtr(9),0,0,0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d3=(betadtr(9)^2-betadtr(3)^2)/sqrt(dcov);
ts_d3;

g=[0,0,0,-2*betadtr(4),0,0,0,0,0,2*betadtr(10),0,0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d4=(betadtr(10)^2-betadtr(4)^2)/sqrt(dcov);
ts_d4;

g=[0,0,0,0,-2*betadtr(5),0,0,0,0,0,2*betadtr(11),0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d5=(betadtr(11)^2-betadtr(5)^2)/sqrt(dcov);
ts_d5;

g=[0,0,0,0,0,-2*betadtr(6),0,0,0,0,0,2*betadtr(12),0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d6=(betadtr(12)^2-betadtr(6)^2)/sqrt(dcov);
ts_d6;
% PLT tr
% Test statistics Variance coefficients (columns 5) 
g=[-2*betaptr(1),0,0,0,0,0,0,0,0,0,0,2*betaptr(12),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p1=(betaptr(12)^2-betaptr(1)^2)/sqrt(dcov);
ts_p1;

g=[0,0,0,0,0,0,-2*betaptr(7),0,0,0,0,0,0,0,0,0,0,2*betaptr(18),0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p2=(betaptr(18)^2-betaptr(7)^2)/sqrt(dcov);
ts_p2;

g=[0,0,0,0,0,0,0,-2*betaptr(8),0,0,0,0,0,0,0,0,0,0,2*betaptr(19),0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p3=(betaptr(19)^2-betaptr(8)^2)/sqrt(dcov);
ts_p3;

g=[0,0,0,0,0,0,0,0,-2*betaptr(9),0,0,0,0,0,0,0,0,0,0,2*betaptr(20),0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p4=(betaptr(20)^2-betaptr(9)^2)/sqrt(dcov);
ts_p4;

g=[0,0,0,0,0,0,0,0,0,-2*betaptr(10),0,0,0,0,0,0,0,0,0,0,2*betaptr(21),0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p5=(betaptr(21)^2-betaptr(10)^2)/sqrt(dcov);
ts_p5;

g=[0,0,0,0,0,0,0,0,0,0,-2*betaptr(11),0,0,0,0,0,0,0,0,0,0,2*betaptr(22),0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p6=(betaptr(22)^2-betaptr(11)^2)/sqrt(dcov);
ts_p6;
% PLT tr
% Test statistics Market impact coefficients (columns 8) 
g=[0,-2*betaptr(2),0,0,0,0,0,0,0,0,0,0,2*betaptr(13),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p21=(betaptr(13)^2-betaptr(2)^2)/sqrt(dcov);
ts_p21;

g=[0,0,-2*betaptr(3),0,0,0,0,0,0,0,0,0,0,2*betaptr(14),0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p31=(betaptr(14)^2-betaptr(3)^2)/sqrt(dcov);
ts_p31;

g=[0,0,0,-2*betaptr(4),0,0,0,0,0,0,0,0,0,0,2*betaptr(15),0,0,0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p41=(betaptr(15)^2-betaptr(4)^2)/sqrt(dcov);
ts_p41;

g=[0,0,0,0,-2*betaptr(5),0,0,0,0,0,0,0,0,0,0,2*betaptr(16),0,0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p51=(betaptr(16)^2-betaptr(5)^2)/sqrt(dcov);
ts_p51;

g=[0,0,0,0,0,-2*betaptr(6),0,0,0,0,0,0,0,0,0,0,2*betaptr(17),0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p61=(betaptr(17)^2-betaptr(6)^2)/sqrt(dcov);
ts_p61;
% PLT tr
% Test statistics Add. Variance coefficients (column 11) 
g=[0,-2*betaptr(7),0,0,0,0,-2*betaptr(2),0,0,0,0,0,2*betaptr(18),0,0,0,0,2*betaptr(13),0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p212=(2*betaptr(13)*betaptr(18)-2*betaptr(2)*betaptr(7))/sqrt(dcov);
ts_p212;

g=[0,0,-2*betaptr(8),0,0,0,0,-2*betaptr(3),0,0,0,0,0,2*betaptr(19),0,0,0,0,2*betaptr(14),0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p313=(2*betaptr(14)*betaptr(19)-2*betaptr(8)*betaptr(3))/sqrt(dcov);
ts_p313;

g=[0,0,0,-2*betaptr(9),0,0,0,0,-2*betaptr(4),0,0,0,0,0,2*betaptr(20),0,0,0,0,2*betaptr(15),0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p414=(2*betaptr(15)*betaptr(20)-2*betaptr(4)*betaptr(9))/sqrt(dcov);
ts_p414;

g=[0,0,0,0,-2*betaptr(10),0,0,0,0,-2*betaptr(5),0,0,0,0,0,2*betaptr(21),0,0,0,0,2*betaptr(16),0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p515=(2*betaptr(16)*betaptr(21)-2*betaptr(5)*betaptr(10))/sqrt(dcov);
ts_p515;

g=[0,0,0,0,0,-2*betaptr(11),0,0,0,0,-2*betaptr(6),0,0,0,0,0,2*betaptr(22),0,0,0,0,2*betaptr(17),0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p616=(2*betaptr(17)*betaptr(22)-2*betaptr(6)*betaptr(11))/sqrt(dcov);
ts_p616;

%% Table 23: Leverage effects in covariance equations
% Scalar tr
% Test statistics Covariance coefficients (column 5) 
g=[-2*betatr(1),2*betatr(2),0]';
dcov=g'*vctr*g;
ts_s=(betatr(2)^2-betatr(1)^2)/sqrt(dcov);
ts_s;
% Diagonal tr
% Test statistics Covariance coefficients (column 5) 
g=[-betadtr(2),-betadtr(1),0,0,0,0,betadtr(8),betadtr(7),0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d12=(betadtr(7)*betadtr(8)-betadtr(1)*betadtr(2))/sqrt(dcov);
ts_d12;

g=[-betadtr(3),0,-betadtr(1),0,0,0,betadtr(9),0,betadtr(7),0,0,0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d13=(betadtr(7)*betadtr(9)-betadtr(1)*betadtr(3))/sqrt(dcov);
ts_d13;

g=[-betadtr(4),0,0,-betadtr(1),0,0,betadtr(10),0,0,betadtr(7),0,0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d14=(betadtr(7)*betadtr(10)-betadtr(1)*betadtr(4))/sqrt(dcov);
ts_d14;

g=[-betadtr(5),0,0,0,-betadtr(1),0,betadtr(11),0,0,0,betadtr(7),0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d15=(betadtr(7)*betadtr(11)-betadtr(1)*betadtr(5))/sqrt(dcov);
ts_d15;

g=[-betadtr(6),0,0,0,0,-betadtr(1),betadtr(12),0,0,0,0,betadtr(7),0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d16=(betadtr(7)*betadtr(12)-betadtr(1)*betadtr(6))/sqrt(dcov);
ts_d16;

g=[0,-betadtr(3),-betadtr(2),0,0,0,0,betadtr(9),betadtr(8),0,0,0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d23=(betadtr(8)*betadtr(9)-betadtr(2)*betadtr(3))/sqrt(dcov);
ts_d23;

g=[0,-betadtr(4),0,-betadtr(2),0,0,0,betadtr(10),0,betadtr(8),0,0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d24=(betadtr(8)*betadtr(10)-betadtr(2)*betadtr(4))/sqrt(dcov);
ts_d24;

g=[0,-betadtr(5),0,0,-betadtr(2),0,0,betadtr(11),0,0,betadtr(8),0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d25=(betadtr(8)*betadtr(11)-betadtr(2)*betadtr(5))/sqrt(dcov);
ts_d25;

g=[0,-betadtr(6),0,0,0,-betadtr(2),0,betadtr(12),0,0,0,betadtr(8),0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d26=(betadtr(8)*betadtr(12)-betadtr(2)*betadtr(6))/sqrt(dcov);
ts_d26;

g=[0,0,-betadtr(4),-betadtr(3),0,0,0,0,betadtr(10),betadtr(9),0,0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d34=(betadtr(9)*betadtr(10)-betadtr(3)*betadtr(4))/sqrt(dcov);
ts_d34;

g=[0,0,-betadtr(5),0,-betadtr(3),0,0,0,betadtr(11),0,betadtr(9),0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d35=(betadtr(9)*betadtr(11)-betadtr(3)*betadtr(5))/sqrt(dcov);
ts_d35;

g=[0,0,-betadtr(6),0,0,-betadtr(3),0,0,betadtr(12),0,0,betadtr(9),0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d36=(betadtr(9)*betadtr(12)-betadtr(3)*betadtr(6))/sqrt(dcov);
ts_d36;

g=[0,0,0,-betadtr(5),-betadtr(4),0,0,0,0,betadtr(11),betadtr(10),0,0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d45=(betadtr(10)*betadtr(11)-betadtr(4)*betadtr(5))/sqrt(dcov);
ts_d45;

g=[0,0,0,-betadtr(6),0,-betadtr(4),0,0,0,betadtr(12),0,betadtr(10),0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d46=(betadtr(10)*betadtr(12)-betadtr(4)*betadtr(6))/sqrt(dcov);
ts_d46;

g=[0,0,0,0,-betadtr(6),-betadtr(5),0,0,0,0,betadtr(12),betadtr(11),0,0,0,0,0,0]';
dcov=g'*vcdtr*g;
ts_d56=(betadtr(11)*betadtr(12)-betadtr(5)*betadtr(6))/sqrt(dcov);
ts_d56;
% PLT tr
% Test statistics Covariance coefficients (column 5) 
g=[-betaptr(7),0,0,0,0,0,-betaptr(1),0,0,0,0,betaptr(18),0,0,0,0,0,betaptr(12),0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p12=(betaptr(12)*betaptr(18)-betaptr(1)*betaptr(7))/sqrt(dcov);
ts_p12;

g=[-betaptr(8),0,0,0,0,0,0,-betaptr(1),0,0,0,betaptr(19),0,0,0,0,0,0,betaptr(12),0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p13=(betaptr(12)*betaptr(19)-betaptr(1)*betaptr(8))/sqrt(dcov);
ts_p13;

g=[-betaptr(9),0,0,0,0,0,0,0,-betaptr(1),0,0,betaptr(20),0,0,0,0,0,0,0,betaptr(12),0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p14=(betaptr(12)*betaptr(20)-betaptr(1)*betaptr(9))/sqrt(dcov);
ts_p14;

g=[-betaptr(10),0,0,0,0,0,0,0,0,-betaptr(1),0,betaptr(21),0,0,0,0,0,0,0,0,betaptr(12),0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p15=(betaptr(12)*betaptr(21)-betaptr(1)*betaptr(10))/sqrt(dcov);
ts_p15;

g=[-betaptr(11),0,0,0,0,0,0,0,0,0,-betaptr(1),betaptr(22),0,0,0,0,0,0,0,0,0,betaptr(12),0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p16=(betaptr(12)*betaptr(22)-betaptr(1)*betaptr(11))/sqrt(dcov);
ts_p16;

g=[0,0,0,0,0,0,-betaptr(8),-betaptr(7),0,0,0,0,0,0,0,0,0,betaptr(19),betaptr(18),0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p23=(betaptr(18)*betaptr(19)-betaptr(7)*betaptr(8))/sqrt(dcov);
ts_p23;

g=[0,0,0,0,0,0,-betaptr(9),0,-betaptr(7),0,0,0,0,0,0,0,0,betaptr(20),0,betaptr(18),0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p24=(betaptr(18)*betaptr(20)-betaptr(7)*betaptr(9))/sqrt(dcov);
ts_p24;

g=[0,0,0,0,0,0,-betaptr(10),0,0,-betaptr(7),0,0,0,0,0,0,0,betaptr(21),0,0,betaptr(18),0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p25=(betaptr(18)*betaptr(21)-betaptr(7)*betaptr(10))/sqrt(dcov);
ts_p25;

g=[0,0,0,0,0,0,-betaptr(11),0,0,0,-betaptr(7),0,0,0,0,0,0,betaptr(22),0,0,0,betaptr(18),0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p26=(betaptr(18)*betaptr(22)-betaptr(7)*betaptr(11))/sqrt(dcov);
ts_p26;

g=[0,0,0,0,0,0,0,-betaptr(9),-betaptr(8),0,0,0,0,0,0,0,0,0,betaptr(20),betaptr(19),0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p34=(betaptr(19)*betaptr(20)-betaptr(8)*betaptr(9))/sqrt(dcov);
ts_p34;

g=[0,0,0,0,0,0,0,-betaptr(10),0,-betaptr(8),0,0,0,0,0,0,0,0,betaptr(21),0,betaptr(19),0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_35=(betaptr(19)*betaptr(21)-betaptr(8)*betaptr(10))/sqrt(dcov);
ts_35;

g=[0,0,0,0,0,0,0,-betaptr(11),0,0,-betaptr(8),0,0,0,0,0,0,0,betaptr(22),0,0,betaptr(19),0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_36=(betaptr(19)*betaptr(22)-betaptr(8)*betaptr(11))/sqrt(dcov);
ts_36;

g=[0,0,0,0,0,0,0,0,-betaptr(10),-betaptr(9),0,0,0,0,0,0,0,0,0,betaptr(21),betaptr(20),0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p45=(betaptr(20)*betaptr(21)-betaptr(9)*betaptr(10))/sqrt(dcov);
ts_p45;

g=[0,0,0,0,0,0,0,0,-betaptr(11),0,-betaptr(9),0,0,0,0,0,0,0,0,betaptr(22),0,betaptr(20),0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p46=(betaptr(20)*betaptr(22)-betaptr(9)*betaptr(11))/sqrt(dcov);
ts_p46;

g=[0,0,0,0,0,0,0,0,0,-betaptr(11),-betaptr(10),0,0,0,0,0,0,0,0,0,betaptr(22),betaptr(21),0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p56=(betaptr(21)*betaptr(22)-betaptr(10)*betaptr(11))/sqrt(dcov);
ts_p56;
% Test statistics Add. Covariance coefficients (columns 8)
g=[-betaptr(2),-betaptr(1),0,0,0,0,0,0,0,0,0,betaptr(13),betaptr(12),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_p121=(betaptr(12)*betaptr(13)-betaptr(1)*betaptr(2))/sqrt(dcov);
ts_p121;

g=[-betaptr(3),0,-betaptr(1),0,0,0,0,0,0,0,0,betaptr(14),0,betaptr(12),0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_131=(betaptr(12)*betaptr(14)-betaptr(1)*betaptr(3))/sqrt(dcov);
ts_131;

g=[-betaptr(4),0,0,-betaptr(1),0,0,0,0,0,0,0,betaptr(15),0,0,betaptr(12),0,0,0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_141=(betaptr(12)*betaptr(15)-betaptr(1)*betaptr(4))/sqrt(dcov);
ts_141;

g=[-betaptr(5),0,0,0,-betaptr(1),0,0,0,0,0,0,betaptr(16),0,0,0,betaptr(12),0,0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_151=(betaptr(12)*betaptr(16)-betaptr(1)*betaptr(5))/sqrt(dcov);
ts_151;

g=[-betaptr(6),0,0,0,0,-betaptr(1),0,0,0,0,0,betaptr(17),0,0,0,0,betaptr(12),0,0,0,0,0,0,0,0,0,0,0]';
dcov=g'*vcptr*g;
ts_161=(betaptr(12)*betaptr(17)-betaptr(1)*betaptr(6))/sqrt(dcov);
ts_161;