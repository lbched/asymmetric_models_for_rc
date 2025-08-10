function [beta,stderr,vc,logl]=Max_lik(lik_fct,b0,vc_type,A,b,Aeq,beq,lb,ub,nonlcon,options,varargin)
 %Computes the maximum-likelihood estimates and associated standard errors
 %-------------------------------------------------------------------------  
 %[beta,stderr,vc,logl]=Max_lik(lik_fct,b0,vc_type,A,b,Aeq,beq,lb,ub,nonlcon,options,varargin);
 %-------------------------------------------------------------------------
 %INPUT: lik_fct - the likelihood function
 %       b0 - vector of initial values of parameters
 %       vc_type - variance-covariance matrix to be implemented
 %                 if vc_type='Sandwich' White method will be used, else numerical Hessian method
 %      
 %       The following inputs are dedicated to the fmincon function:
 %       A - Vector, constrains for the minimization
 %       b - Vector, constrains for the minimization
 %       Aeq - Vector, constrains for the minimization
 %       beq - Vector, constrains for the minimization
 %       lb - Vector, lower bound for the parameters
 %       ub - Vector, upper bound for the parameters
 %       nonlcon - subjects the minimization to the nonlinear inequalities
 %       options - function's options
 %       varargin - input arguments that are needed to calculate the likelihood function 
 %-------------------------------------------------------------------------
 %OUTPUT: beta - Vector of estimated parameters
 %        stderr - Vector of associated standard errors
 %        vc - Variance covariance matrix
 %        logl - Value of the log-likelihood
 %-------------------------------------------------------------------------
f0=feval(lik_fct,b0,varargin{:});
T=size(f0,1);

if T==1
    error('Likelihood function should return a column vector of scores')
end

[beta,fval,~,~,~,~,~] =...
    fmincon(@ml_sum,b0,A,b,Aeq,beq,lb,ub,nonlcon,options,lik_fct,varargin{:});

hessian=HessMp(@ml_sum,beta,lik_fct,varargin{:});

inv_h=-hessian\eye(size(hessian,1));

if strcmp(vc_type,'Sandwich')
    
% disp('estimation of the variance-covariance matrix via Sandwich [White]');
g  = gradp(lik_fct,beta,varargin{:});
vc = inv_h*(g'*g)*inv_h;
    
else % default is information matrix
    
% disp('estimation of the variance-covariance matrix via Hessian'); 
vc = -inv_h;
    
end

stderr  = sqrt(diag(vc));
logl    = -fval;

function l=ml_sum(b,lik_fct,varargin)
l=feval(lik_fct,b,varargin{:});
l=-sum(l);

function g=gradp(f,x0,varargin) % computes the gradient of f evaluated at x
% f should return either a scalar or a column vector
% x0 should be a column vector of parameters
f0=feval(f,x0,varargin{:}); 
[T,~]=size(f0);

if size(x0,2)>size(x0,1)
x0=x0';
end
k=size(x0,1); % number of parameters wrt which one should compute the gradient

h=0.0000001; % some small number

g=zeros(T,k); 
e=eye(k); 
for j=1:k
    if x0(j)>1 % if argument is big enough, compute relative number   
        f1=feval(f,(x0.*( ones(k,1) +  e(:,j) *h )),varargin{:});    
        g(:,j)=(f1-f0)/(x0(j)*h);    
    else
        f1=feval(f, x0 +  e(:,j) *h ,varargin{:});    
        g(:,j)=(f1-f0)/h;    
    
    end
    
end

function H=HessMp(f,x0,varargin) % computes the Hessian matrix of f evaluated at x0
% f should return either a scalar or a column vector
% x0 should be a column vector of parameters
f0=feval(f,x0,varargin{:}); 
[~,co]=size(f0);
if co>1; error('Error in HessMp, The function should be a column vector or a scalar'); end

[k,c]=size(x0);
if k<c
x0=x0';
end
k=size(x0,1); % number of parameters wrt which one should compute the gradient

h=0.00001; % some small number

H=zeros(k,k); % will contain the Hessian
e=eye(k); 

h2=h/2;
for ii=1:k
    if x0(ii)>100 % if argument is big enough, compute relative number   
        x0P= x0.*( ones(k,1) +  e(:,ii) *h2 );
        x0N= x0.*( ones(k,1) -  e(:,ii) *h2 );
        Deltaii = x0(ii)*h;
    else
        x0P = x0 +  e(:,ii) *h2;
        x0N = x0 -  e(:,ii) *h2;
        Deltaii = h;
    end
    
    for jj=1:ii
    if x0(jj)>100 % if argument is big enough, compute relative number   
        x0PP = x0P .* ( ones(k,1) +  e(:,jj) *h2 );
        x0PN = x0P .* ( ones(k,1) -  e(:,jj) *h2 );
        x0NP = x0N .* ( ones(k,1) +  e(:,jj) *h2 );
        x0NN = x0N .* ( ones(k,1) -  e(:,jj) *h2 );
        Delta = Deltaii*x0(jj)*h;
    else
        x0PP = x0P  +  e(:,jj) *h2; 
        x0PN = x0P  -  e(:,jj) *h2; 
        x0NP = x0N  +  e(:,jj) *h2; 
        x0NN = x0N  -  e(:,jj) *h2; 
        Delta = Deltaii*h;
    end
    
        fPP = feval(f,x0PP,varargin{:});  % forward,forward
        fPN = feval(f,x0PN,varargin{:});  % forward,backward
        fNP = feval(f,x0NP,varargin{:});  % backward,forward
        fNN = feval(f,x0NN,varargin{:});  % backward,backward
        
        H(ii,jj)=(sum(fPP)-sum(fPN)-sum(fNP)+sum(fNN))/Delta;
        H(jj,ii)=H(ii,jj);
    end
end