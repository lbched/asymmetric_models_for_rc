%% Forecasting Results (Tables:10-11)
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
%% Selected Forecast Generations
%% Selected Scalar models
% sym
M=76;
T1=2137;
betasym=zeros(2,380);
stderrsym=zeros(2,380);
vcsym=zeros(2,2,380);
loglsym=zeros(380,1);
A=zeros(380,1); 
B=zeros(380,1);
Qt=zeros(n,n,380);
for j=1:5
lb = zeros(2,1)';
ub = ones(2,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Scsymllk(p,RC(:,:,1+j*M:T1+j*M),n,T1);
x0 = [0.1; 0.8]';
[betasym(:,j),stderrsym(:,j),vcsym(:,:,j),loglsym(j)] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
    A(j)=betasym(1,j);
    B(j)=betasym(2,j);
Qt(:,:,1:M)=Scsymrec(A(1),B(1),RC(:,:,T1+1:T1+M-1),n,M); 
Qt(:,:,M+1:2*M)=Scsymrec(A(2),B(2),RC(:,:,T1+M:T1+2*M-1),n,M);
Qt(:,:,2*M+1:3*M)=Scsymrec(A(3),B(3),RC(:,:,T1+2*M:T1+3*M-1),n,M);
Qt(:,:,3*M+1:4*M)=Scsymrec(A(4),B(4),RC(:,:,T1+3*M:T1+4*M-1),n,M);
Qt(:,:,4*M+1:5*M)=Scsymrec(A(5),B(5),RC(:,:,T1+4*M:T1+5*M-1),n,M);
end
symFor=zeros(380,n*(n+1)/2);
for i=1:380
symFor(i,:)=vech(Qt(:,:,i),n);
end
% tr
M=76;
T1=2137;
betatr=zeros(3,380);
stderrtr=zeros(3,380);
vctr=zeros(3,3,380);
logltr=zeros(380,1);
A1=zeros(380,1); 
A2=zeros(380,1);
B=zeros(380,1);
Qt=zeros(n,n,380);
for j=1:5
lb = zeros(3,1)';
ub = ones(3,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Sctrllk(p,RC(:,:,1+j*M:T1+j*M),PMC(:,:,1+j*M:T1+j*M),NC(:,:,1+j*M:T1+j*M),n,T1);
x0 = [0.1; 0.07; 0.8]';
[betatr(:,j),stderrtr(:,j),vctr(:,:,j),logltr(j)] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
    A1(j)=betatr(1,j);
    A2(j)=betatr(2,j);
    B(j)=betatr(3,j);
Qt(:,:,1:M)=Sctrrec(A1(1),A2(1),B(1),RC(:,:,T1+1:T1+M-1),PMC(:,:,T1+1:T1+M-1),NC(:,:,T1+1:T1+M-1),n,M); 
Qt(:,:,M+1:2*M)=Sctrrec(A1(2),A2(2),B(2),RC(:,:,T1+M:T1+2*M-1),PMC(:,:,T1+M:T1+2*M-1),NC(:,:,T1+M:T1+2*M-1),n,M);
Qt(:,:,2*M+1:3*M)=Sctrrec(A1(3),A2(3),B(3),RC(:,:,T1+2*M:T1+3*M-1),PMC(:,:,T1+2*M:T1+3*M-1),NC(:,:,T1+2*M:T1+3*M-1),n,M);
Qt(:,:,3*M+1:4*M)=Sctrrec(A1(4),A2(4),B(4),RC(:,:,T1+3*M:T1+4*M-1),PMC(:,:,T1+3*M:T1+4*M-1),NC(:,:,T1+3*M:T1+4*M-1),n,M);
Qt(:,:,4*M+1:5*M)=Sctrrec(A1(5),A2(5),B(5),RC(:,:,T1+4*M:T1+5*M-1),PMC(:,:,T1+4*M:T1+5*M-1),NC(:,:,T1+4*M:T1+5*M-1),n,M);
end
trFor=zeros(380,n*(n+1)/2);
for i=1:380
trFor(i,:)=vech(Qt(:,:,i),n);
end
% tr^oc
M=76;
T1=2137;
betatr=zeros(3,380);
stderrtr=zeros(3,380);
vctr=zeros(3,3,380);
logltr=zeros(380,1);
A1=zeros(380,1); 
A2=zeros(380,1);
B=zeros(380,1);
Qt=zeros(n,n,380);
for j=1:5
lb = zeros(3,1)';
ub = ones(3,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Sctrllk(p,RC(:,:,1+j*M:T1+j*M),PMC_oc(:,:,1+j*M:T1+j*M),NC_oc(:,:,1+j*M:T1+j*M),n,T1);
x0 = [0.1; 0.07; 0.8]';
[betatr(:,j),stderrtr(:,j),vctr(:,:,j),logltr(j)] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
    A1(j)=betatr(1,j);
    A2(j)=betatr(2,j);
    B(j)=betatr(3,j);
Qt(:,:,1:M)=Sctrrec(A1(1),A2(1),B(1),RC(:,:,T1+1:T1+M-1),PMC_oc(:,:,T1+1:T1+M-1),NC_oc(:,:,T1+1:T1+M-1),n,M); 
Qt(:,:,M+1:2*M)=Sctrrec(A1(2),A2(2),B(2),RC(:,:,T1+M:T1+2*M-1),PMC_oc(:,:,T1+M:T1+2*M-1),NC_oc(:,:,T1+M:T1+2*M-1),n,M);
Qt(:,:,2*M+1:3*M)=Sctrrec(A1(3),A2(3),B(3),RC(:,:,T1+2*M:T1+3*M-1),PMC_oc(:,:,T1+2*M:T1+3*M-1),NC_oc(:,:,T1+2*M:T1+3*M-1),n,M);
Qt(:,:,3*M+1:4*M)=Sctrrec(A1(4),A2(4),B(4),RC(:,:,T1+3*M:T1+4*M-1),PMC_oc(:,:,T1+3*M:T1+4*M-1),NC_oc(:,:,T1+3*M:T1+4*M-1),n,M);
Qt(:,:,4*M+1:5*M)=Sctrrec(A1(5),A2(5),B(5),RC(:,:,T1+4*M:T1+5*M-1),PMC_oc(:,:,T1+4*M:T1+5*M-1),NC_oc(:,:,T1+4*M:T1+5*M-1),n,M);
end
tr_ocFor=zeros(380,n*(n+1)/2);
for i=1:380
tr_ocFor(i,:)=vech(Qt(:,:,i),n);
end
%% Selected Diagonal models
% sym
n=6;
M=76;
T1=2137;
betasym=zeros(12,5);
stderrsym=zeros(12,5);
vcsym=zeros(12,12,5);
loglsym=zeros(5,1);
A=zeros(6,6,5); 
B=zeros(6,6,5);
Qt=zeros(6,6,380);
for j=1:5
lb=zeros(n*2,1)';
ub=ones(n*2,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Diagsymllk(p,RC(:,:,1+j*M:2137+j*M),n,T1);
x0=[sqrt(0.27)*ones(n,1); sqrt(0.7)*ones(n,1)]';
[betasym(:,j),stderrsym(:,j),vcsym(:,:,j),loglsym(j)] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
A(:,:,j)=diag((betasym(1:n,j)));
B(:,:,j)=diag((betasym((n+1):(2*n),j)));
Q1 = vech(mean(RC(:,:,1:2517),3),n);
Qt(:,:,1:M)=Diagsymrec1(A(:,:,1),B(:,:,1),RC(:,:,T1+1:T1+M-1),n,M,Q1); 
Qt(:,:,M+1:2*M)=Diagsymrec1(A(:,:,2),B(:,:,2),RC(:,:,T1+M:T1+2*M-1),n,M,Q1);
Qt(:,:,2*M+1:3*M)=Diagsymrec1(A(:,:,3),B(:,:,3),RC(:,:,T1+2*M:T1+3*M-1),n,M,Q1);
Qt(:,:,3*M+1:4*M)=Diagsymrec1(A(:,:,4),B(:,:,4),RC(:,:,T1+3*M:T1+4*M-1),n,M,Q1);
Qt(:,:,4*M+1:5*M)=Diagsymrec1(A(:,:,5),B(:,:,5),RC(:,:,T1+4*M:T1+5*M-1),n,M,Q1);
end
dsymFor=zeros(380,n*(n+1)/2);
for i=1:380
dsymFor(i,:)=vech(Qt(:,:,i),n);
end
% tr
M=76;
T1=2137;
betatr=zeros(18,5);
stderrtr=zeros(18,5);
vctr=zeros(18,18,5);
logltr=zeros(5,1);
A1=zeros(6,6,5); 
A2=zeros(6,6,5);
B=zeros(6,6,5);
Qt=zeros(6,6,380);
for j=1:5
lb = zeros(18,1)';
ub = ones(18,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Diagtrllk(p,RC(:,:,1+j*M:2137+j*M),PMC(:,:,1+j*M:2137+j*M),NC(:,:,1+j*M:2137+j*M),n,T1);
x0=[sqrt(0.12)*ones(n,1); sqrt(0.14)*ones(n,1); sqrt(0.8)*ones(n,1)]';
[betatr(:,j),stderrtr(:,j),vctr(:,:,j),logltr(j)] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
A1(:,:,j)=diag((betatr(1:n,j)));
A2(:,:,j)=diag((betatr((n+1):(2*n),j)));
B(:,:,j)=diag((betatr((2*n+1):(3*n),j)));
Q1 = vech(mean(RC(:,:,1:2517),3),n);
P=kron(mean(PMC(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2),mean(PMC(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2));
N1=kron(mean(NC(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2),mean(NC(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2));
Qt(:,:,1:M)=Diagtrrec1(A1(:,:,1),A2(:,:,1),B(:,:,1),PMC(:,:,T1+1:T1+M-1),NC(:,:,T1+1:T1+M-1),n,M,Q1,P,N1); 
Qt(:,:,M+1:2*M)=Diagtrrec1(A1(:,:,2),A2(:,:,2),B(:,:,2),PMC(:,:,T1+M:T1+2*M-1),NC(:,:,T1+M:T1+2*M-1),n,M,Q1,P,N1); 
Qt(:,:,2*M+1:3*M)=Diagtrrec1(A1(:,:,3),A2(:,:,3),B(:,:,3),PMC(:,:,T1+2*M:T1+3*M-1),NC(:,:,T1+2*M:T1+3*M-1),n,M,Q1,P,N1); 
Qt(:,:,3*M+1:4*M)=Diagtrrec1(A1(:,:,4),A2(:,:,4),B(:,:,4),PMC(:,:,T1+3*M:T1+4*M-1),NC(:,:,T1+3*M:T1+4*M-1),n,M,Q1,P,N1); 
Qt(:,:,4*M+1:5*M)=Diagtrrec1(A1(:,:,5),A2(:,:,5),B(:,:,5),PMC(:,:,T1+4*M:T1+5*M-1),NC(:,:,T1+4*M:T1+5*M-1),n,M,Q1,P,N1); 
end
dtrFor=zeros(380,21);
for i=1:380
dtrFor(i,:)=vech(Qt(:,:,i),n);
end
% tr^oc
M=76;
betatr=zeros(18,5);
stderrtr=zeros(18,5);
vctr=zeros(18,18,5);
logltr=zeros(5,1);
A1=zeros(6,6,5); 
A2=zeros(6,6,5);
B=zeros(6,6,5);
Qt=zeros(6,6,380);
for j=1:5
lb = zeros(18,1)';
ub = ones(18,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) Diagtrllk(p,RC(:,:,1+j*M:2137+j*M),PMC_oc(:,:,1+j*M:2137+j*M),NC_oc(:,:,1+j*M:2137+j*M),n,T1);
x0=[sqrt(0.12)*ones(n,1); sqrt(0.14)*ones(n,1); sqrt(0.8)*ones(n,1)]';
[betatr(:,j),stderrtr(:,j),vctr(:,:,j),logltr(j)] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
A1(:,:,j)=diag((betatr(1:n,j)));
A2(:,:,j)=diag((betatr((n+1):(2*n),j)));
B(:,:,j)=diag((betatr((2*n+1):(3*n),j)));
Q1 = vech(mean(RC(:,:,1:2517),3),n);
P=kron(mean(PMC_oc(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2),mean(PMC_oc(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2));
N1=kron(mean(NC_oc(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2),mean(NC_oc(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2));
Qt(:,:,1:M)=Diagtrrec1(A1(:,:,1),A2(:,:,1),B(:,:,1),PMC_oc(:,:,T1+1:T1+M-1),NC_oc(:,:,T1+1:T1+M-1),n,M,Q1,P,N1); 
Qt(:,:,M+1:2*M)=Diagtrrec1(A1(:,:,2),A2(:,:,2),B(:,:,2),PMC_oc(:,:,T1+M:T1+2*M-1),NC_oc(:,:,T1+M:T1+2*M-1),n,M,Q1,P,N1); 
Qt(:,:,2*M+1:3*M)=Diagtrrec1(A1(:,:,3),A2(:,:,3),B(:,:,3),PMC_oc(:,:,T1+2*M:T1+3*M-1),NC_oc(:,:,T1+2*M:T1+3*M-1),n,M,Q1,P,N1); 
Qt(:,:,3*M+1:4*M)=Diagtrrec1(A1(:,:,4),A2(:,:,4),B(:,:,4),PMC_oc(:,:,T1+3*M:T1+4*M-1),NC_oc(:,:,T1+3*M:T1+4*M-1),n,M,Q1,P,N1); 
Qt(:,:,4*M+1:5*M)=Diagtrrec1(A1(:,:,5),A2(:,:,5),B(:,:,5),PMC_oc(:,:,T1+4*M:T1+5*M-1),NC_oc(:,:,T1+4*M:T1+5*M-1),n,M,Q1,P,N1); 
end
dtr_ocFor=zeros(380,21);
for i=1:380
dtr_ocFor(i,:)=vech(Qt(:,:,i),n);
end
%% Selected PLT models
% sym
M=76;
betasym=zeros(17,5);
stderrsym=zeros(17,5);
vcsym=zeros(17,17,5);
loglsym=zeros(5,1);
A=zeros(6,6,5); 
B=zeros(6,6,5);
Qt=zeros(6,6,380);
for j=1:5
lb=[0 -0.3 -0.3 -0.3 -0.3 -0.3 0 0 0 0 0 0 0 0 0 0 0]';
ub=ones(17,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) PLTsymllk(p,RC(:,:,1+j*M:2137+j*M),n,T1);
x0=[0.13;0.02;0.02;0.02;0.02;0.02;0.15;0.15;0.15;0.15;0.15;0.5;0.5;0.5;0.5;0.5;0.5]';
[betasym(:,j),stderrsym(:,j),vcsym(:,:,j),loglsym(j)] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
A(:,:,j)=eye(6);
A(1,1,j)=(betasym(1,j));
A(2,1,j)=(betasym(2,j));
A(3,1,j)=(betasym(3,j));
A(4,1,j)=(betasym(4,j));
A(5,1,j)=(betasym(5,j));
A(6,1,j)=(betasym(6,j));
A(2,2,j)=(betasym(7,j));
A(3,3,j)=(betasym(8,j));
A(4,4,j)=(betasym(9,j));
A(5,5,j)=(betasym(10,j));
A(6,6,j)=(betasym(11,j));
B(:,:,j)=eye(6);
B(1,1,j)=(betasym(12,j));
B(2,2,j)=(betasym(13,j));
B(3,3,j)=(betasym(14,j));
B(4,4,j)=(betasym(15,j));
B(5,5,j)=(betasym(16,j));
B(6,6,j)=(betasym(17,j));
Q1 = vech(mean(RC(:,:,1:2517),3),n);
Qt(:,:,1:M)=PLTsymrec1(A(:,:,1),B(:,:,1),RC(:,:,T1+1:T1+M-1),n,M,Q1); 
Qt(:,:,M+1:2*M)=PLTsymrec1(A(:,:,2),B(:,:,2),RC(:,:,T1+M:T1+2*M-1),n,M,Q1);
Qt(:,:,2*M+1:3*M)=PLTsymrec1(A(:,:,3),B(:,:,3),RC(:,:,T1+2*M:T1+3*M-1),n,M,Q1);
Qt(:,:,3*M+1:4*M)=PLTsymrec1(A(:,:,4),B(:,:,4),RC(:,:,T1+3*M:T1+4*M-1),n,M,Q1);
Qt(:,:,4*M+1:5*M)=PLTsymrec1(A(:,:,5),B(:,:,5),RC(:,:,T1+4*M:T1+5*M-1),n,M,Q1);
end
psymFor=zeros(380,21);
for i=1:380
psymFor(i,:)=vech(Qt(:,:,i),n);
end
% tr
M=76;
betatr=zeros(28,5);
stderrtr=zeros(28,5);
vctr=zeros(28,28,5);
logltr=zeros(5,1);
A1=zeros(6,6,5); 
A2=zeros(6,6,5);
B=zeros(6,6,5);
Qt=zeros(6,6,380);
for j=1:5
lb=[0 -0.1 -0.1 -0.1 -0.1 -0.1 0 0 0 0 0 0 -0.1 -0.1 -0.1 -0.1 -0.1 0 0 0 0 0 0 0 0 0 0 0]';
ub=ones(28,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) PLTtrllk(p,RC(:,:,1+j*M:2137+j*M),PMC(:,:,1+j*M:2137+j*M),NC(:,:,1+j*M:2137+j*M),n,T1);
x0=[0.368804558362593;0.01;0.01;0.01;0.01;0.01;	0.538065952083568;	0.530318973022856;	0.511549377512276;	0.545485168848220;	0.573257572781041;	0.472004744259667;0.01;0.01;0.01;0.01;0.01;	0.567506080656019;	0.554722139862587;	0.533979918990755;	0.577588629369931;	0.611362715636051;	0.896864395332714;	0.791443957839432;	0.809195896886628;	0.815725618142746;	0.786526806403777;	0.766984572255703]';
[betatr(:,j),stderrtr(:,j),vctr(:,:,j),logltr(j)] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
A1(:,:,j)=eye(6);
A1(1,1,j)=(betatr(1,j));
A1(2,1,j)=(betatr(2,j));
A1(3,1,j)=(betatr(3,j));
A1(4,1,j)=(betatr(4,j));
A1(5,1,j)=(betatr(5,j));
A1(6,1,j)=(betatr(6,j));
A1(2,2,j)=(betatr(7,j));
A1(3,3,j)=(betatr(8,j));
A1(4,4,j)=(betatr(9,j));
A1(5,5,j)=(betatr(10,j));
A1(6,6,j)=(betatr(11,j));
A2(:,:,j)=eye(6);
A2(1,1,j)=(betatr(12,j));
A2(2,1,j)=(betatr(13,j));
A2(3,1,j)=(betatr(14,j));
A2(4,1,j)=(betatr(15,j));
A2(5,1,j)=(betatr(16,j));
A2(6,1,j)=(betatr(17,j));
A2(2,2,j)=(betatr(18,j));
A2(3,3,j)=(betatr(19,j));
A2(4,4,j)=(betatr(20,j));
A2(5,5,j)=(betatr(21,j));
A2(6,6,j)=(betatr(22,j));
B(:,:,j)=eye(6);
B(1,1,j)=(betatr(23,j));
B(2,2,j)=(betatr(24,j));
B(3,3,j)=(betatr(25,j));
B(4,4,j)=(betatr(26,j));
B(5,5,j)=(betatr(27,j));
B(6,6,j)=(betatr(28,j));
Q1 = vech(mean(RC(:,:,1:2517),3),n);
P=kron(mean(PMC(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2),mean(PMC(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2));
N1=kron(mean(NC(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2),mean(NC(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2));
Qt(:,:,1:M)=PLTtrrec1(A1(:,:,1),A2(:,:,1),B(:,:,1),PMC(:,:,T1+1:T1+M-1),NC(:,:,T1+1:T1+M-1),n,M,Q1,P,N1); 
Qt(:,:,M+1:2*M)=PLTtrrec1(A1(:,:,2),A2(:,:,2),B(:,:,2),PMC(:,:,T1+M:T1+2*M-1),NC(:,:,T1+M:T1+2*M-1),n,M,Q1,P,N1); 
Qt(:,:,2*M+1:3*M)=PLTtrrec1(A1(:,:,3),A2(:,:,3),B(:,:,3),PMC(:,:,T1+2*M:T1+3*M-1),NC(:,:,T1+2*M:T1+3*M-1),n,M,Q1,P,N1); 
Qt(:,:,3*M+1:4*M)=PLTtrrec1(A1(:,:,4),A2(:,:,4),B(:,:,4),PMC(:,:,T1+3*M:T1+4*M-1),NC(:,:,T1+3*M:T1+4*M-1),n,M,Q1,P,N1); 
Qt(:,:,4*M+1:5*M)=PLTtrrec1(A1(:,:,5),A2(:,:,5),B(:,:,5),PMC(:,:,T1+4*M:T1+5*M-1),NC(:,:,T1+4*M:T1+5*M-1),n,M,Q1,P,N1); 
end
ptrFor=zeros(380,21);
for i=1:380
ptrFor(i,:)=vech(Qt(:,:,i),n);
end
% tr^oc
M=76;
betatr=zeros(28,5);
stderrtr=zeros(28,5);
vctr=zeros(28,28,5);
logltr=zeros(5,1);
A1=zeros(6,6,5); 
A2=zeros(6,6,5);
B=zeros(6,6,5);
Qt=zeros(6,6,380);
for j=1:5
lb=[0 -0.1 -0.1 -0.1 -0.1 -0.1 0 0 0 0 0 0 -0.1 -0.1 -0.1 -0.1 -0.1 0 0 0 0 0 0 0 0 0 0 0]';
ub=ones(28,1)';
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
llk = @(p) PLTtrllk(p,RC(:,:,1+j*M:2137+j*M),PMC_oc(:,:,1+j*M:2137+j*M),NC_oc(:,:,1+j*M:2137+j*M),n,T1);
x0=[0.368804558362593;0.01;0.01;0.01;0.01;0.01;	0.538065952083568;	0.530318973022856;	0.511549377512276;	0.545485168848220;	0.573257572781041;	0.472004744259667;0.01;0.01;0.01;0.01;0.01;	0.567506080656019;	0.554722139862587;	0.533979918990755;	0.577588629369931;	0.611362715636051;	0.896864395332714;	0.791443957839432;	0.809195896886628;	0.815725618142746;	0.786526806403777;	0.766984572255703]';
[betatr(:,j),stderrtr(:,j),vctr(:,:,j),logltr(j)] = Max_lik(llk,x0,'Sandwich',[],[],[],[],lb,ub,[],options);
A1(:,:,j)=eye(6);
A1(1,1,j)=(betatr(1,j));
A1(2,1,j)=(betatr(2,j));
A1(3,1,j)=(betatr(3,j));
A1(4,1,j)=(betatr(4,j));
A1(5,1,j)=(betatr(5,j));
A1(6,1,j)=(betatr(6,j));
A1(2,2,j)=(betatr(7,j));
A1(3,3,j)=(betatr(8,j));
A1(4,4,j)=(betatr(9,j));
A1(5,5,j)=(betatr(10,j));
A1(6,6,j)=(betatr(11,j));
A2(:,:,j)=eye(6);
A2(1,1,j)=(betatr(12,j));
A2(2,1,j)=(betatr(13,j));
A2(3,1,j)=(betatr(14,j));
A2(4,1,j)=(betatr(15,j));
A2(5,1,j)=(betatr(16,j));
A2(6,1,j)=(betatr(17,j));
A2(2,2,j)=(betatr(18,j));
A2(3,3,j)=(betatr(19,j));
A2(4,4,j)=(betatr(20,j));
A2(5,5,j)=(betatr(21,j));
A2(6,6,j)=(betatr(22,j));
B(:,:,j)=eye(6);
B(1,1,j)=(betatr(23,j));
B(2,2,j)=(betatr(24,j));
B(3,3,j)=(betatr(25,j));
B(4,4,j)=(betatr(26,j));
B(5,5,j)=(betatr(27,j));
B(6,6,j)=(betatr(28,j));
Q1 = vech(mean(RC(:,:,1:2517),3),n);
P=kron(mean(PMC_oc(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2),mean(PMC_oc(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2));
N1=kron(mean(NC_oc(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2),mean(NC_oc(:,:,1:2517),3)^(1/2)*mean(RC(:,:,1:2517),3)^(-1/2));
Qt(:,:,1:M)=PLTtrrec1(A1(:,:,1),A2(:,:,1),B(:,:,1),PMC_oc(:,:,T1+1:T1+M-1),NC_oc(:,:,T1+1:T1+M-1),n,M,Q1,P,N1); 
Qt(:,:,M+1:2*M)=PLTtrrec1(A1(:,:,2),A2(:,:,2),B(:,:,2),PMC_oc(:,:,T1+M:T1+2*M-1),NC_oc(:,:,T1+M:T1+2*M-1),n,M,Q1,P,N1); 
Qt(:,:,2*M+1:3*M)=PLTtrrec1(A1(:,:,3),A2(:,:,3),B(:,:,3),PMC_oc(:,:,T1+2*M:T1+3*M-1),NC_oc(:,:,T1+2*M:T1+3*M-1),n,M,Q1,P,N1); 
Qt(:,:,3*M+1:4*M)=PLTtrrec1(A1(:,:,4),A2(:,:,4),B(:,:,4),PMC_oc(:,:,T1+3*M:T1+4*M-1),NC_oc(:,:,T1+3*M:T1+4*M-1),n,M,Q1,P,N1); 
Qt(:,:,4*M+1:5*M)=PLTtrrec1(A1(:,:,5),A2(:,:,5),B(:,:,5),PMC_oc(:,:,T1+4*M:T1+5*M-1),NC_oc(:,:,T1+4*M:T1+5*M-1),n,M,Q1,P,N1); 
end
ptr_ocFor=zeros(380,21);
for i=1:380
ptr_ocFor(i,:)=vech(Qt(:,:,i),n);
end

%% The file paths are relative to the folder Forecasts/
m = 9; % total number of models per version
T = 2517; % sample size
T1 = 2137; % estimation window
T2 = 380; % out-of-sample window

% Scalar forecasts
% sym
A1 = symFor;
F1s = zeros(n,n,T2);
for t = 1:T2
   F1s(:,:,t) = buildSymmetric(A1(t,:),n);
end
% tr
A2 = trFor;
F2s = zeros(n,n,T2);
for t = 1:T2
   F2s(:,:,t) = buildSymmetric(A2(t,:),n);
end
% trPNM 
A3 = readmatrix("strPNM.csv");
F3s = zeros(n,n,T2);
for t = 1:T2
   F3s(:,:,t) = buildSymmetric(A3(t,:),n);
end
% trPNtauM
A4 = readmatrix("strPNtauM.csv");
F4s = zeros(n,n,T2);
for t = 1:T2
   F4s(:,:,t) = buildSymmetric(A4(t,:),n);
end
% semi
A5 = readmatrix("ssemi.csv");
F5s = zeros(n,n,T2);
for t = 1:T2
   F5s(:,:,t) = buildSymmetric(A5(t,:),n);
end
% semi-tau
A6 = readmatrix("ssemitau.csv");
F6s = zeros(n,n,T2);
for t = 1:T2
   F6s(:,:,t) = buildSymmetric(A6(t,:),n);
end
% tr_oc
A7 = tr_ocFor;
F7s = zeros(n,n,T2);
for t = 1:T2
   F7s(:,:,t) = buildSymmetric(A7(t,:),n);
end
% trPNM_oc
A8 = readmatrix("strPNM^oc.csv");
F8s = zeros(n,n,T2);
for t = 1:T2
   F8s(:,:,t) = buildSymmetric(A8(t,:),n);
end 
% trPNtauM_oc
A9 = readmatrix("strPNtauM^oc.csv");
F9s = zeros(n,n,T2);
for t = 1:T2
   F9s(:,:,t) = buildSymmetric(A9(t,:),n);
end

% Diagonal forecasts
% sym
A1 = dsymFor;
F1d = zeros(n,n,T2);
for t = 1:T2
   F1d(:,:,t) = buildSymmetric(A1(t,:),n);
end
% tr
A2 = dtrFor;
F2d = zeros(n,n,T2);
for t = 1:T2
   F2d(:,:,t) = buildSymmetric(A2(t,:),n);
end
% trPNM
A3 = readmatrix("dtrPNM.csv");
F3d = zeros(n,n,T2);
for t = 1:T2
   F3d(:,:,t) = buildSymmetric(A3(t,:),n);
end
% trPNtauM
A4 = readmatrix("dtrPNtauM.csv");
F4d = zeros(n,n,T2);
for t = 1:T2
   F4d(:,:,t) = buildSymmetric(A4(t,:),n);
end
% semi
A5 = readmatrix("dsemi.csv");
F5d = zeros(n,n,T2);
for t = 1:T2
   F5d(:,:,t) = buildSymmetric(A5(t,:),n);
end
% semi-tau
A6 = readmatrix("dsemitau.csv");
F6d = zeros(n,n,T2);
for t = 1:T2
   F6d(:,:,t) = buildSymmetric(A6(t,:),n);
end
% tr_oc
A7 = dtr_ocFor;
F7d = zeros(n,n,T2);
for t = 1:T2
   F7d(:,:,t) = buildSymmetric(A7(t,:),n);
end
% trPNM_oc
A8 = readmatrix("dtrPNM^oc.csv");
F8d = zeros(n,n,T2);
for t = 1:T2
   F8d(:,:,t) = buildSymmetric(A8(t,:),n);
end 
% trPNtauM_oc
A9 = readmatrix("dtrPNtauM^oc.csv");
F9d = zeros(n,n,T2);
for t = 1:T2
   F9d(:,:,t) = buildSymmetric(A9(t,:),n);
end

% PLT forecasts
% sym
A1 = psymFor;
F1p = zeros(n,n,T2);
for t = 1:T2
   F1p(:,:,t) = buildSymmetric(A1(t,:),n);
end
% tr
A2 = ptrFor;
F2p = zeros(n,n,T2);
for t = 1:T2
   F2p(:,:,t) = buildSymmetric(A2(t,:),n);
end
% trPNM
A3 = readmatrix("plttrPNM.csv");
F3p = zeros(n,n,T2);
for t = 1:T2
   F3p(:,:,t) = buildSymmetric(A3(t,:),n);
end
% trPNtauM
A4 = readmatrix("plttrPNtauM.csv");
F4p = zeros(n,n,T2);
for t = 1:T2
   F4p(:,:,t) = buildSymmetric(A4(t,:),n);
end
% semi
A5 = readmatrix("pltsemi.csv");
F5p = zeros(n,n,T2);
for t = 1:T2
   F5p(:,:,t) = buildSymmetric(A5(t,:),n);
end
% semi-tau
A6 = readmatrix("pltsemitau.csv");
F6p = zeros(n,n,T2);
for t = 1:T2
   F6p(:,:,t) = buildSymmetric(A6(t,:),n);
end
% tr_oc
A7 = ptr_ocFor;
F7p = zeros(n,n,T2);
for t = 1:T2
   F7p(:,:,t) = buildSymmetric(A7(t,:),n);
end
% trPNM_oc
A8 = readmatrix("plttrPNM^oc.csv");
F8p = zeros(n,n,T2);
for t = 1:T2
   F8p(:,:,t) = buildSymmetric(A8(t,:),n);
end 
% trPNtauM_oc
A9 = readmatrix("plttrPNtauM^oc.csv");
F9p = zeros(n,n,T2);
for t = 1:T2
   F9p(:,:,t) = buildSymmetric(A9(t,:),n);
end

% RC as an Unobserved Covariance proxy to compute loss functions
rc = readmatrix('RC.csv');
p = zeros(n,n,T2);
for t = 1:T2
   p(:,:,t) = buildSymmetric(rc((t+T1),:),n)*25200;
end

%% The file paths are relative to the folder Fun/.
%% Table 10: FN loss (columns 5-6)
% Scalar models
LF_FNs = zeros(T2,m);      % matrix for the loss functions
for t=1:T2
    LF_FNs(t,1)=sqrt(trace((p(:,:,t)-F1s(:,:,t))'*(p(:,:,t)-F1s(:,:,t))));  % sym
    LF_FNs(t,2)=sqrt(trace((p(:,:,t)-F2s(:,:,t))'*(p(:,:,t)-F2s(:,:,t))));  % tr
    LF_FNs(t,3)=sqrt(trace((p(:,:,t)-F3s(:,:,t))'*(p(:,:,t)-F3s(:,:,t))));  % trPNM
    LF_FNs(t,4)=sqrt(trace((p(:,:,t)-F4s(:,:,t))'*(p(:,:,t)-F4s(:,:,t))));  % trPNtauM
    LF_FNs(t,5)=sqrt(trace((p(:,:,t)-F7s(:,:,t))'*(p(:,:,t)-F7s(:,:,t))));  % tr_oc
    LF_FNs(t,6)=sqrt(trace((p(:,:,t)-F8s(:,:,t))'*(p(:,:,t)-F8s(:,:,t))));  % trPNM_oc
    LF_FNs(t,7)=sqrt(trace((p(:,:,t)-F9s(:,:,t))'*(p(:,:,t)-F9s(:,:,t))));  % trPNtauM_oc
    LF_FNs(t,8)=sqrt(trace((p(:,:,t)-F5s(:,:,t))'*(p(:,:,t)-F5s(:,:,t))));  % semi
    LF_FNs(t,9)=sqrt(trace((p(:,:,t)-F6s(:,:,t))'*(p(:,:,t)-F6s(:,:,t))));  % semi-tau
end
% Diagonal models
LF_FNd = zeros(T2,m);      % matrix for the loss functions
for t=1:T2
    LF_FNd(t,1)=sqrt(trace((p(:,:,t)-F1d(:,:,t))'*(p(:,:,t)-F1d(:,:,t))));  % sym
    LF_FNd(t,2)=sqrt(trace((p(:,:,t)-F2d(:,:,t))'*(p(:,:,t)-F2d(:,:,t))));  % tr
    LF_FNd(t,3)=sqrt(trace((p(:,:,t)-F3d(:,:,t))'*(p(:,:,t)-F3d(:,:,t))));  % trPNM
    LF_FNd(t,4)=sqrt(trace((p(:,:,t)-F4d(:,:,t))'*(p(:,:,t)-F4d(:,:,t))));  % trPNtauM
    LF_FNd(t,5)=sqrt(trace((p(:,:,t)-F7d(:,:,t))'*(p(:,:,t)-F7d(:,:,t))));  % tr_oc
    LF_FNd(t,6)=sqrt(trace((p(:,:,t)-F8d(:,:,t))'*(p(:,:,t)-F8d(:,:,t))));  % trPNM_oc
    LF_FNd(t,7)=sqrt(trace((p(:,:,t)-F9d(:,:,t))'*(p(:,:,t)-F9d(:,:,t))));  % trPNtauM_oc
    LF_FNd(t,8)=sqrt(trace((p(:,:,t)-F5d(:,:,t))'*(p(:,:,t)-F5d(:,:,t))));  % semi
    LF_FNd(t,9)=sqrt(trace((p(:,:,t)-F6d(:,:,t))'*(p(:,:,t)-F6d(:,:,t))));  % semi-tau
end
% PLT models
LF_FNp = zeros(T2,m);      % matrix for the loss functions
for t=1:T2
    LF_FNp (t,1)=sqrt(trace((p(:,:,t)-F1p(:,:,t))'*(p(:,:,t)-F1p(:,:,t))));  % sym
    LF_FNp (t,2)=sqrt(trace((p(:,:,t)-F2p(:,:,t))'*(p(:,:,t)-F2p(:,:,t))));  % tr
    LF_FNp (t,3)=sqrt(trace((p(:,:,t)-F3p(:,:,t))'*(p(:,:,t)-F3p(:,:,t))));  % trPNM
    LF_FNp (t,4)=sqrt(trace((p(:,:,t)-F4p(:,:,t))'*(p(:,:,t)-F4p(:,:,t))));  % trPNtauM
    LF_FNp (t,5)=sqrt(trace((p(:,:,t)-F7p(:,:,t))'*(p(:,:,t)-F7p(:,:,t))));  % tr_oc
    LF_FNp (t,6)=sqrt(trace((p(:,:,t)-F8p(:,:,t))'*(p(:,:,t)-F8p(:,:,t))));  % trPNM_oc
    LF_FNp (t,7)=sqrt(trace((p(:,:,t)-F9p(:,:,t))'*(p(:,:,t)-F9p(:,:,t))));  % trPNtauM_oc
    LF_FNp (t,8)=sqrt(trace((p(:,:,t)-F5p(:,:,t))'*(p(:,:,t)-F5p(:,:,t))));  % semi
    LF_FNp (t,9)=sqrt(trace((p(:,:,t)-F6p(:,:,t))'*(p(:,:,t)-F6p(:,:,t))));  % semi-tau
end
%% Average FN scalar models (column 5)
AVFNs=zeros(1,m);
for j=1:m
    AVFNs(:,j)=mean(LF_FNs(:,j));
end
AVFNs;
%% Average FN diagonal models (column 5)
AVFNd=zeros(1,m);
for j=1:m
    AVFNd(:,j)=mean(LF_FNd(:,j));
end
AVFNd;
%% Average FN PLT models (column 5)
AVFNp=zeros(1,m);
for j=1:m
   AVFNp(:,j)=mean(LF_FNp(:,j));
end
AVFNp;
% To compute the resulting MCSs, adopt the functions mcs.m & statiomary_bootstrap.m from the MFE Toolbox of Kevin Sheppard.
LF_FN=[LF_FNs LF_FNd LF_FNp];
[includedFN,~,~,~,~,~]=mcs(LF_FN,0.10,10000,100); % MCS90% 
includedFN;
% includedFN: 20(PLT tr)

%% Table 11: GMVP loss unconstrained (columns 3-4)
% Scalar models
% sym
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma1s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F1s(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F1s(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma1s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma2s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F2s(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F2s(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma2s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma3s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F3s(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F3s(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma3s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma4s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F4s(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F4s(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma4s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma5s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F5s(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F5s(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma5s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi-tau
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma6s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F6s(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F6s(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma6s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma7s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F7s(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F7s(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma7s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma8s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F8s(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F8s(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma8s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma9s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F9s(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F9s(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma9s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% Diagonal models
% sym
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma1d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F1d(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F1d(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma1d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma2d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F2d(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F2d(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma2d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma3d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F3d(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F3d(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma3d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma4d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F4d(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F4d(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma4d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma5d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F5d(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F5d(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma5d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi-tau
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma6d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F6d(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F6d(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma6d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma7d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F7d(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F7d(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma7d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma8d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F8d(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F8d(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma8d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma9d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F9d(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F9d(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma9d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% PLT models
% sym
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma1p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F1p(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F1p(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma1p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma2p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F2p(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F2p(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma2p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma3p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F3p(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F3p(:,:,t),n),6)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma3p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma4p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F4p(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F4p(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma4p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma5p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F5p(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F5p(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma5p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi-tau
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma6p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F6p(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F6p(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma6p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma7p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F7p(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F7p(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma7p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma8p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F8p(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F8p(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma8p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma9p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F9p(:,:,t),n),n),'AssetMean',zeros(1,n));
prob = optimproblem('ObjectiveSense','minimize');
x = optimvar('x',n,1,'LowerBound',-1); 
prob.Objective = x'*buildSymmetric(vech(F9p(:,:,t),n),n)*x;
prob.Constraints.sumToTau = sum(x) == 1;
sol = solve(prob);
w(t,:) = sol.x;
sigma9p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
%% Average GMVP scalar models (column 3)
LF_GMVs = zeros(T2,m);      % matrix for the loss functions 
for t=1:T2
    LF_GMVs(t,1) = sigma1s(t);  % sym
    LF_GMVs(t,2) = sigma2s(t);  % tr
    LF_GMVs(t,3) = sigma3s(t);  % trPNM
    LF_GMVs(t,4) = sigma4s(t);  % trPNtauM
    LF_GMVs(t,5) = sigma7s(t);  % tr_oc
    LF_GMVs(t,6) = sigma8s(t);  % trPNM_oc
    LF_GMVs(t,7) = sigma9s(t);  % trPNtauM_oc
    LF_GMVs(t,8) = sigma5s(t);  % semi
    LF_GMVs(t,9) = sigma6s(t);  % semi-tau
end
AVGMVs = zeros(1,9);
for j = 1:9
    AVGMVs(:,j) = mean(LF_GMVs(:,j));
end
AVGMVs;
%% Average GMVP diagonal models (column 3)
LF_GMVd = zeros(T2,m);      % matrix for the loss functions 
for t=1:T2
    LF_GMVd(t,1) = sigma1d(t);  % sym
    LF_GMVd(t,2) = sigma2d(t);  % tr
    LF_GMVd(t,3) = sigma3d(t);  % trPNM
    LF_GMVd(t,4) = sigma4d(t);  % trPNtauM
    LF_GMVd(t,5) = sigma7d(t);  % tr_oc
    LF_GMVd(t,6) = sigma8d(t);  % trPNM_oc
    LF_GMVd(t,7) = sigma9d(t);  % trPNtauM_oc
    LF_GMVd(t,8) = sigma5d(t);  % semi
    LF_GMVd(t,9) = sigma6d(t);  % semi-tau
end
AVGMVd = zeros(1,9);
for j = 1:9
    AVGMVd(:,j) = mean(LF_GMVd(:,j));
end
AVGMVd;
%% Average GMVP PLT models (column 3)
LF_GMVp = zeros(T2,m);      % matrix for the loss functions 
for t=1:T2
    LF_GMVp(t,1) = sigma1p(t);  % sym
    LF_GMVp(t,2) = sigma2p(t);  % tr
    LF_GMVp(t,3) = sigma3p(t);  % trPNM
    LF_GMVp(t,4) = sigma4p(t);  % trPNtauM
    LF_GMVp(t,5) = sigma7p(t);  % tr_oc
    LF_GMVp(t,6) = sigma8p(t);  % trPNM_oc
    LF_GMVp(t,7) = sigma9p(t);  % trPNtauM_oc
    LF_GMVp(t,8) = sigma5p(t);  % semi
    LF_GMVp(t,9) = sigma6p(t);  % semi-tau
end
AVGMVp = zeros(1,9);
for j = 1:9
    AVGMVp(:,j) = mean(LF_GMVp(:,j));
end
AVGMVp;
% To compute the resulting MCSs, adopt the functions mcs.m & statiomary_bootstrap.m from the MFE Toolbox of Kevin Sheppard.
LF_GMV=[LF_GMVs LF_GMVd LF_GMVp];
[includedGMV,~,~,~,~,~]=mcs(LF_GMV,0.10,10000,100); % MCS90% 
includedGMV;
% includedGMV: 2(scalar tr),5(scalar tr_oc),1(scalar sym),7(scalar trPNtauM_oc),6(scalar trPNM_oc),
%              3(scalar trPNM),4(scalar trPNtauM)

%% Table 11: GMVP loss constrained (columns 5-6)
% Scalar models
% sym
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma1cs = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F1s(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma1cs(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma2cs = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F2s(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma2cs(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma3cs = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F3s(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma3cs(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma4cs = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F4s(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma4cs(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma5cs = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F5s(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma5cs(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi-tau
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma6cs = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F6s(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma6cs(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma7cs = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F7s(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma7cs(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma8cs = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F8s(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma8cs(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma9cs = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F9s(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma9cs(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% Diagonal models
% sym
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma1cd = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F1d(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma1cd(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma2cd = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F2d(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma2cd(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma3cd = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F3d(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma3cd(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma4cd = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F4d(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma4cd(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma5cd = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F5d(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma5cd(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi-tau
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma6cd = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F6d(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma6cd(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma7cd = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F7d(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma7cd(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma8cd = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F8d(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma8cd(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma9cd = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F9d(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma9cd(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% PLT models
% sym
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma1cp = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F1p(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma1cp(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma2cp = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F2p(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma2cp(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma3cp = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F3p(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma3cp(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma4cp = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F4p(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma4cp(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma5cp = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F5p(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma5cp(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi-tau
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma6cp = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F6p(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma6cp(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma7cp = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F7p(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma7cp(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma8cp = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F8p(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma8cp(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM_oc
pt = repmat(Portfolio,T2,1);
w = zeros(T2,n);
sigma9cp = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F9p(:,:,t),n),n),'AssetMean',zeros(1,n));
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierLimits(pt(t,:),'min');
sigma9cp(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
%% Average constrained GMVP scalar models (column 5)
LFGVs = zeros(T2,m);      % matrix for the loss functions 
for t=1:T2
    LFGVs(t,1) = sigma1cs(t);  % sym
    LFGVs(t,2) = sigma2cs(t);  % tr
    LFGVs(t,3) = sigma3cs(t);  % trPNM
    LFGVs(t,4) = sigma4cs(t);  % trPNtauM
    LFGVs(t,5) = sigma7cs(t);  % tr_oc
    LFGVs(t,6) = sigma8cs(t);  % trPNM_oc
    LFGVs(t,7) = sigma9cs(t);  % trPNtauM_oc
    LFGVs(t,8) = sigma5cs(t);  % semi
    LFGVs(t,9) = sigma6cs(t);  % semi-tau
end
AVGMVcs = zeros(1,m);
for j = 1:m
    AVGMVcs(:,j) = mean(LFGVs(:,j));
end
AVGMVcs;
%% Average constrained GMVP diagonal models (column 5)
LFGVd = zeros(T2,m);      % matrix for the loss functions 
for t=1:T2
    LFGVd(t,1) = sigma1cd(t);  % sym
    LFGVd(t,2) = sigma2cd(t);  % tr
    LFGVd(t,3) = sigma3cd(t);  % trPNM
    LFGVd(t,4) = sigma4cd(t);  % trPNtauM
    LFGVd(t,5) = sigma7cd(t);  % tr_oc
    LFGVd(t,6) = sigma8cd(t);  % trPNM_oc
    LFGVd(t,7) = sigma9cd(t);  % trPNtauM_oc
    LFGVd(t,8) = sigma5cd(t);  % semi
    LFGVd(t,9) = sigma6cd(t);  % semi-tau
end
AVGMVcd = zeros(1,m);
for j = 1:m
    AVGMVcd(:,j) = mean(LFGVd(:,j));
end
AVGMVcd;
%% Average constrained GMVP PLT models (column 5)
LFGVp = zeros(T2,m);      % matrix for the loss functions 
for t=1:T2
    LFGVp(t,1) = sigma1cp(t);  % sym
    LFGVp(t,2) = sigma2cp(t);  % tr
    LFGVp(t,3) = sigma3cp(t);  % trPNM
    LFGVp(t,4) = sigma4cp(t);  % trPNtauM
    LFGVp(t,5) = sigma7cp(t);  % tr_oc
    LFGVp(t,6) = sigma8cp(t);  % trPNM_oc
    LFGVp(t,7) = sigma9cp(t);  % trPNtauM_oc
    LFGVp(t,8) = sigma5cp(t);  % semi
    LFGVp(t,9) = sigma6cp(t);  % semi-tau
end
AVGMVcp = zeros(1,m);
for j = 1:m
    AVGMVcp(:,j) = mean(LFGVp(:,j));
end
AVGMVcp;
% To compute the resulting MCSs, adopt the functions mcs.m & statiomary_bootstrap.m from the MFE Toolbox of Kevin Sheppard.
LFGMVc=[LFGVs LFGVd LFGVp];
[includedGMVc,~,~,~,~,~]=mcs(LFGMVc,0.10,10000,100); % MCS90% 
includedGMVc;
% includedGMVc: 2(scalar tr)

%% Table 11: MVP loss constrained (columns 7-8)
% Scalar models
% sym
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma11s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F1s(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma11s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma21s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F2s(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma21s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma31s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F3s(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma31s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma41s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F4s(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma41s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma51s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F5s(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma51s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi-tau
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma61s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F6s(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma61s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr_oc
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma71s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F7s(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma71s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM_oc
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma81s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F8s(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma81s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM_oc
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma91s = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F9s(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma91s(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% Diagonal models
% sym
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma11d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F1d(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma11d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma21d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F2d(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma21d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma31d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F3d(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma31d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma41d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F4d(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma41d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma51d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F5d(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma51d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi-tau
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma61d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F6d(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma61d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr_oc
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma71d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F7d(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma71d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM_oc
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma81d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F8d(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma81d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM_oc
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma91d = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F9d(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma91d(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% PLT models
% sym
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma11p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F1p(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma11p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma21p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F2p(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma21p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma31p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F3p(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma31p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma41p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F4p(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma41p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma51p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F5p(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma51p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% semi-tau
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma61p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F6p(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma61p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% tr_oc
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma71p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F7p(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma71p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNM_oc
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma81p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F8p(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma81p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
% trPNtauM_oc
pt = repmat(Portfolio,T1,1);
w = zeros(T2,n);
sigma91p = zeros(T2,1);
for t = 1:T2
pt(t,:) = Portfolio('AssetCovar',buildSymmetric(vech(F9p(:,:,t),n),n),'AssetMean',[0.00000001, 0.000002, 0.000001, 0.000002, 0.000001, 0.000002]');
pt(t,:) = setDefaultConstraints(pt(t,:));
w(t,:) = estimateFrontierByReturn(pt(t,:),0.0000015);
sigma91p(t)=sqrt(w(t,:)*p(:,:,t)*w(t,:)');
end
%% Average MVP scalar models (column 7)
LFMVs = zeros(T2,m);      % matrix for the loss functions 
for t=1:T2
    LFMVs(t,1) = sigma11s(t);  % sym
    LFMVs(t,2) = sigma21s(t);  % tr
    LFMVs(t,3) = sigma31s(t);  % trPNM
    LFMVs(t,4) = sigma41s(t);  % trPNtauM
    LFMVs(t,5) = sigma71s(t);  % tr_oc
    LFMVs(t,6) = sigma81s(t);  % trPNM_oc
    LFMVs(t,7) = sigma91s(t);  % trPNtauM_oc
    LFMVs(t,8) = sigma51s(t);  % semi
    LFMVs(t,9) = sigma61s(t);  % semi-tau
end
AVMVs = zeros(1,9);
for j = 1:9
    AVMVs(:,j) = mean(LFMVs(:,j));
end
AVMVs;
%% Average MVP diagonal models (column 7)
LFMVd = zeros(T2,m);      % matrix for the loss functions 
for t=1:T2
    LFMVd(t,1) = sigma11d(t);  % sym
    LFMVd(t,2) = sigma21d(t);  % tr
    LFMVd(t,3) = sigma31d(t);  % trPNM
    LFMVd(t,4) = sigma41d(t);  % trPNtauM
    LFMVd(t,5) = sigma71d(t);  % tr_oc
    LFMVd(t,6) = sigma81d(t);  % trPNM_oc
    LFMVd(t,7) = sigma91d(t);  % trPNtauM_oc
    LFMVd(t,8) = sigma51d(t);  % semi
    LFMVd(t,9) = sigma61d(t);  % semi-tau
end
AVMVd = zeros(1,m);
for j = 1:m
    AVMVd(:,j) = mean(LFMVd(:,j));
end
AVMVd;
%% Average MVP PLT models (column 7)
LFMVp = zeros(T2,m);      % matrix for the loss functions 
for t=1:T2
    LFMVp(t,1) = sigma11p(t);  % sym
    LFMVp(t,2) = sigma21p(t);  % tr
    LFMVp(t,3) = sigma31p(t);  % trPNM
    LFMVp(t,4) = sigma41p(t);  % trPNtauM
    LFMVp(t,5) = sigma71p(t);  % tr_oc
    LFMVp(t,6) = sigma81p(t);  % trPNM_oc
    LFMVp(t,7) = sigma91p(t);  % trPNtauM_oc
    LFMVp(t,8) = sigma51p(t);  % semi
    LFMVp(t,9) = sigma61p(t);  % semi-tau
end
AVMVp = zeros(1,m);
for j = 1:m
    AVMVp(:,j) = mean(LFMVp(:,j));
end
AVMVp;
% To compute the resulting MCSs, adopt the functions mcs.m & statiomary_bootstrap.m from the MFE Toolbox of Kevin Sheppard.
LFMV=[LFMVs LFMVd LFMVp];
[includedMV,~,~,~,~,~]=mcs(LFMV,0.10,10000,50); % MCS90% 
includedMV;
% includedMV: 26(PLT semi),7(scalar trPNtauM_oc),2(scalar tr),5(scalar tr_oc),6(scalar trPNM_oc),
%             1(scalar sym),8(scalar semi),25(PLT trPNtauM_oc),9(scalar semi-tau)