function [U]=PCexpansion_onevariable(marg_mean,marg_var,Distributiontype)

% clear;clc;
% mean=100;variance=400;order=7;Distributiontype=1;

VarianceV=[1 1 2 6 24 120 720 5040];
%VarianceV=[1 1 2 6 24 120 720 5040 40320 362880];
range=8;

switch Distributiontype
    case 1                                        %Gamma Distribution
        
        beta=marg_var/marg_mean; alpha=marg_mean/beta;      %shape and scale parameter
        U(1)=integral(@(z) icdf('Gamma',cdf('Normal',z,0,1),alpha,beta).*(1).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(1),-range,range);
        U(2)=integral(@(z) icdf('Gamma',cdf('Normal',z,0,1),alpha,beta).*(z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(2),-range,range);
        U(3)=integral(@(z) icdf('Gamma',cdf('Normal',z,0,1),alpha,beta).*(z.^2-1).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(3),-range,range);
        U(4)=integral(@(z) icdf('Gamma',cdf('Normal',z,0,1),alpha,beta).*(z.^3-3*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(4),-range,range);
        U(5)=integral(@(z) icdf('Gamma',cdf('Normal',z,0,1),alpha,beta).*(z.^4-6*z.^2+3).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(5),-range,range);
        U(6)=integral(@(z) icdf('Gamma',cdf('Normal',z,0,1),alpha,beta).*(z.^5-10*z.^3+15*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(6),-range,range);
        U(7)=integral(@(z) icdf('Gamma',cdf('Normal',z,0,1),alpha,beta).*(z.^6-15*z.^4+45*z.^2-15).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(7),-range,range);
        U(8)=integral(@(z) icdf('Gamma',cdf('Normal',z,0,1),alpha,beta).*(z.^7-21*z.^5+105*z.^3-105*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(8),-range,range);
           
    case 2                                         %lognormal Distribution
        
        sigma=sqrt(log(marg_var/marg_mean^2+1)); miu=log(marg_mean)-sigma^2/2;    %log mean and log standard deviation
        U(1)=integral(@(z) icdf('Lognormal',cdf('Normal',z,0,1),miu,sigma).*(1).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(1),-range,range);
        U(2)=integral(@(z) icdf('Lognormal',cdf('Normal',z,0,1),miu,sigma).*(z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(2),-range,range);
        U(3)=integral(@(z) icdf('Lognormal',cdf('Normal',z,0,1),miu,sigma).*(z.^2-1).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(3),-range,range);
        U(4)=integral(@(z) icdf('Lognormal',cdf('Normal',z,0,1),miu,sigma).*(z.^3-3*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(4),-range,range);
        U(5)=integral(@(z) icdf('Lognormal',cdf('Normal',z,0,1),miu,sigma).*(z.^4-6*z.^2+3).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(5),-range,range);
        U(6)=integral(@(z) icdf('Lognormal',cdf('Normal',z,0,1),miu,sigma).*(z.^5-10*z.^3+15*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(6),-range,range);
        U(7)=integral(@(z) icdf('Lognormal',cdf('Normal',z,0,1),miu,sigma).*(z.^6-15*z.^4+45*z.^2-15).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(7),-range,range);
        U(8)=integral(@(z) icdf('Lognormal',cdf('Normal',z,0,1),miu,sigma).*(z.^7-21*z.^5+105*z.^3-105*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(8),-range,range);
                 
    case 3                                          %Guass Distribution
        miu=marg_mean; sigma=sqrt(marg_var);             %gauss mean and standard deviation  
        U(1)=integral(@(z) icdf('Normal',cdf('Normal',z,0,1),miu,sigma).*(1).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(1),-range,range);
        U(2)=integral(@(z) icdf('Normal',cdf('Normal',z,0,1),miu,sigma).*(z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(2),-range,range);
        U(3)=integral(@(z) icdf('Normal',cdf('Normal',z,0,1),miu,sigma).*(z.^2-1).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(3),-range,range);
        U(4)=integral(@(z) icdf('Normal',cdf('Normal',z,0,1),miu,sigma).*(z.^3-3*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(4),-range,range);
        U(5)=integral(@(z) icdf('Normal',cdf('Normal',z,0,1),miu,sigma).*(z.^4-6*z.^2+3).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(5),-range,range);
        U(6)=integral(@(z) icdf('Normal',cdf('Normal',z,0,1),miu,sigma).*(z.^5-10*z.^3+15*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(6),-range,range);
        U(7)=integral(@(z) icdf('Normal',cdf('Normal',z,0,1),miu,sigma).*(z.^6-15*z.^4+45*z.^2-15).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(7),-range,range);
        U(8)=integral(@(z) icdf('Normal',cdf('Normal',z,0,1),miu,sigma).*(z.^7-21*z.^5+105*z.^3-105*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(8),-range,range);
        
    case 4                         %distribution with user-provided PDF numerical values
        % needs more verification, use with caution
        pdf_x=marg_mean;
        pdf_y=marg_var;
        
        cdf_x=pdf_x;
        cdf_y=cumtrapz(pdf_x,pdf_y);
        cdf_y(1)=cdf_y(2)/2;
        increment=0.0001;
        U=zeros(1,8);
        range_left=icdf('Normal',cdf_y(1),0,1);
        range_right=icdf('Normal',cdf_y(end),0,1);
        
        for z=range_left:increment:range_right
            
            inv_value=interp1(cdf_y,cdf_x, cdf('Normal',z,0,1), 'linear','extrap');

            U(1)=U(1)+inv_value*(1).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(1)*increment;
            U(2)=U(2)+inv_value*(z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(2)*increment;
            U(3)=U(3)+inv_value*(z.^2-1).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(3)*increment;
            U(4)=U(4)+inv_value*(z.^3-3*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(4)*increment;
            U(5)=U(5)+inv_value*(z.^4-6*z.^2+3).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(5)*increment;
            U(6)=U(6)+inv_value*(z.^5-10*z.^3+15*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(6)*increment;
            U(7)=U(7)+inv_value*(z.^6-15*z.^4+45*z.^2-15).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(7)*increment;
            U(8)=U(8)+inv_value*(z.^7-21*z.^5+105*z.^3-105*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(8)*increment;
%           U(9)=U(9)+inv_value*(z.^8-28*z.^6+210*z.^4-420*z.^2+105).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(9)*increment;
%           U(10)=U(10)+inv_value*(z.^9-36*z.^7+378*z.^5-1260*z.^3+945*z).*exp(-z.^2/2)/sqrt(2*pi)/VarianceV(10)*increment;
        end
        
end
%check variance
check_var=sum(VarianceV(2:end).*(U(2:end).*U(2:end)));


%%check PDF of soil shear modulus
% y = pdf('Gamma',0:1:500,2,50); plot(0:1:500,y); hold on;
% y = pdf('Lognormal',0:1:500,miu,sigma); plot(0:1:500,y);
% legend('Gamma distribution','Lognormal distribution');
% title('use two distributions for same soil uncertainty (mean=100, variance=5000)');

