%%
%Authors: Fangbo Wang (fbwang@ucdavis.edu), Hexiang Wang(hexwang@ucdavis.edu)

%Date-Aug21-2018
%Note that this following code can be used to quantify a non-Gaussian
%non-stationary random process or random field using Hermite polynomial chaos.
%The correlation of random process/field is captured by KL expansion.

%The implemented algorithm can be referred to 
%"Sakamoto, S., Ghanem, R. Polynomial Chaos Decomposition 
%for the Simulation of Non-Gaussian Nonstationary Stochastic Processes. 
%Journal of Engineering Mechanics, 2002". 


clear;clc;

%%
%------------------------------------------------------------------------
%              Part 1.1: geometry and mesh of the random field/process
%------------------------------------------------------------------------
totallength=10; leftbound=0; rightbound=leftbound+totallength; %1D model 
NodeNum=11; %total number of KL nodes in the mesh
Le=totallength/(NodeNum-1);  %element length, two-noded
xx=leftbound:Le:rightbound;  %coordinate of nodes

%Marginal information (can be any pdf function with user-uploaded "pdf_x" and "pdf_y")
%0--user-provided PDF data; 1--Gamma dist; 2--Lognormal dist; 3--Gauss dist;  
% load('pdf.mat');    % upload user-defined pdf data
Distributiontype=1;      
%marginal mean of random field/process, it can vary along space or time
marg_mean=100*ones(NodeNum,1); 
%marginal variance of random field/process, it can vary along space or time
marg_var=1600*ones(NodeNum,1);  


%Correlation information (covar_mat can also be any function with user-uploaded "covar_mat")
% load('covar_mat.mat');  % upload user-defined covar_mat
lc=5;                %correlation length
Correlationtype=1;   %0-user defined correlation function, 1-exponential function; 2-square exponential function
[covar_mat]=Construct_covar(xx, NodeNum, marg_var, Correlationtype, lc); % construct the covariance matrix


%%
%-------------------------------------------------------------------------
%    Part 1.2: Polynomial expansion information for the random field/process
%-------------------------------------------------------------------------
% set required Dimension and order for the PC and load them
Dim=10; order_p=4;              %PC dimension and order
load('PC_Variance_dim10order4.mat');
P=Num_PCterms(Dim,order_p);  %number of PC terms 
PC=PC(1:P,:);  Variance=Variance(1:P);


%%
%Part 2.2: One dimensional PC expansion of a random variable
for i=1:NodeNum
    
    if Distributiontype<=3
       [U]=PCexpansion_onevariable(marg_mean(i), marg_var(i), Distributiontype); % for Distributiontype=1,2,3
        
    elseif Distributiontype==4
        [U]=PCexpansion_onevariable(pdf_x, pdf_y, Distributiontype); % for Distributiontype=4
        
    end
    PCweights(i,:)=U;  % U is PC coeffs for 1D PC expansion
end

%Part 2.3: Inversely compute the basis (standard Gaussian) correlation in terms of 
%          the given covariance structure of random field/process per eq. 20,21 

%Comment: inverse order is ideally to be equal to PC order, however, it is
%         simplified to be order 1 to avoid solving the nonliner equation 
%         which could render no reasonable root or more than two reasonable roots.
inverse_order=4;
[covar_gaussian]=Inverse_covar(inverse_order, NodeNum, covar_mat, PCweights);


%%%%%%verification only----synthesized covariance of input random field/process 
% syn_covar=zeros(NodeNum);
% for i=1:NodeNum
%     for j=1:NodeNum
%         for order=1:inverse_order
%             syn_covar(i,j) = syn_covar(i,j) + PCweights(i,order+1)*PCweights(j,order+1)*factorial(order)*covar_gaussian(i,j)^order;
%         end
%     end
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 1D KL-expansion of standard Gaussian random field/process %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%comment: For a 1D KL-expansion, KL is a simple eigenvalue problem.

[V,D] = eig(covar_gaussian);

% since Matlab's eig can't return ordered eigenvalue, sorting is needed.
[D_val,index] = sort(diag(D), 'descend');
D = D(index,index);
V = V(:,index);


%%%% smooth the eigenvector to avoid numerical jumps
%%%% usually it is not needed
% for i=1:NodeNum
%     pp1=csaps(xx(2:end-1),V(2:end-1,i),1);
%     y1=ppval(pp1,xx);
%     V(:,i)=y1;
% end

% Normalization of the eigenvector
for i=1:NodeNum
    Norm=sqrt(sum(V(:,i).^2));
    V(:,i)=V(:,i)/Norm;
end

%%%% verification only---correlation matrix of gaussian process
% syn_cov_gaussian=zeros(NodeNum,NodeNum);
% for i=1:NodeNum
%     for j=1:NodeNum
%         for k=1:NodeNum
%             syn_cov_gaussian(i,j)=syn_cov_gaussian(i,j)+D(k,k)*V(i,k)*V(j,k);
%         end
%     end
% end

% Computation of the normalized factor for each node location
for Location=1:NodeNum
    temp=0;
    for i=1:NodeNum
        temp=temp+D(i,i)*V(Location,i)^2;
    end
    factor(Location)=sqrt(1/temp);
end

%%%%% truncation of KL-expansion, denomi should be put on denominator of PC coeffs computation
denomi=zeros(NodeNum,1);
for i=1:NodeNum
    for j=1:Dim
       denomi(i)=denomi(i)+D(j,j)*(factor(i)*V(i,j))^2;
    end
end
denomi=sqrt(denomi);


%%
%Part 2.5: calculate coefficient of multidimensional expansion at each node per equation 16 of Sakamoto and Ghanem
u=zeros(NodeNum,size(PC,1));
u(:,1)=U(1);
P=size(PC,1); % number of PC terms
for k=1:NodeNum                       %k denote location of node
    for i=2:size(PC,1)

        p=sum(PC(i,2:2:end));
        u(k,i)=factorial(p)*U(p+1)/Variance(i);
        numx=zeros(Dim,1);
        for ii=1:2:size(PC,2)
            if PC(i,ii)~=0    
                numx(PC(i,ii))=PC(i,ii+1);    
            end
        end
        
        % introduction of tmep_u is necessary for efficiency, it avoid the
        % repeated moving of pointer to obtain u(k,i) in the forloop
        temp_u=u(k,i); 
        for j=1:Dim
              temp_u=temp_u*(sqrt(D(j,j))*factor(k)*V(k,j)/denomi(k))^numx(j);
        end
        u(k,i)=temp_u;
    end  
    k
end
PC_coeffs=u;


%%

%check mean and SD with actual
synMean=PC_coeffs(:,1);
synVar=zeros(NodeNum,1);
for i=1:NodeNum
    for j=2:P
        synVar(i)=synVar(i)+PC_coeffs(i,j)^2*Variance(j);
    end
end

%check correlation
correlation=zeros(NodeNum);
for Loc1=1:NodeNum
    for Loc2=1:NodeNum
        for i=2:P
            correlation(Loc1,Loc2)=correlation(Loc1,Loc2)+u(Loc1,i)*u(Loc2,i)*Variance(i)/sqrt(synVar(Loc1))/sqrt(synVar(Loc2));
        end
    end
end

for i=1:NodeNum
    temp=cat(1,diag(correlation,i-1),diag(correlation,1-i));
    syn_corr(i)=mean(temp);
end

%%
%compare mean, variane, correlation
figure
subplot(3,1,1)
plot(synMean,'b'); hold on;
plot(marg_mean,'r'); grid on;
ylim([max(marg_mean-sqrt(marg_var))  max(marg_mean+sqrt(marg_var))]);
xlabel('lag time'); ylabel('mean');
legend('synthesized mean','exact mean recordings');
title('comparison of mean,SD,correaltion with Dim2order1(164 KL elements)');

subplot(3,1,2)
plot(sqrt(synVar),'b'); hold on;
plot(sqrt(marg_var),'r'); grid on;
ylim([0 2*max(sqrt(marg_var))]);
xlabel('lag time'); ylabel('SD');
legend('synthesized SD','exact SD recordings');

subplot(3,1,3)
% surf(correlation); hold on;
% mesh(correlation);
plot(xx,syn_corr,'b--'); hold on;
plot(xx,exp(-xx/lc),'b-'); grid on;
xlabel('LAG DISTANCE'); ylabel('COEF. OF CORRELATION');
legend('Synthesized corr.','Exact coor.');


% corr=correlation(2:5:end,2:5:end);
% tol=abs(corr_acc_327by327(1:326,1:326)-corr);
% mesh(tol);
% legend('error of corr');

