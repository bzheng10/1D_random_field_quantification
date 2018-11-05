
function [covar_mat]=Construct_covar(xx, NodeNum, marg_var, Correlationtype, lc)
    
    covar_mat=zeros(NodeNum);
    
    switch Correlationtype
        case 0  %user-defined correlation function
            
                %please upload the correlation matrix
                load('corr_mat.mat', 'corr_mat');
                for i=1:NodeNum
                    for j=1:NodeNum
                        covar_mat(i,j)=sqrt(marg_var(i))*sqrt(marg_var(j))*corr_mat(i,j);
                    end
                end
                
                
        case 1  %exponential function
            
                for i=1:NodeNum
                    for j=1:NodeNum
                        covar_mat(i,j)=sqrt(marg_var(i))*sqrt(marg_var(j))*exp(-abs(xx(i)-xx(j))/lc);
                    end
                end
                
        case 2  %square exponential function
            
                for i=1:NodeNum
                    for j=1:NodeNum
                        covar_mat(i,j)=sqrt(marg_var(i))*sqrt(marg_var(j))*exp(-abs(xx(i)-xx(j))^2/lc^2);
                    end
                end

    end

end

