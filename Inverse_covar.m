
function [covar_gaussian]=Inverse_covar(inverse_order, NodeNum, covar_mat, PCweights)

   
    covar_gaussian=zeros(NodeNum);
    switch inverse_order
        case 1
            % to simplify to order 1, it is a linear equation
            for i=1:NodeNum
                for j=1:NodeNum
                    covar_gaussian(i,j)=covar_mat(i,j)/PCweights(i,2)/PCweights(j,2);
                end
            end
        otherwise
            % for higher-order, a nonlinear equation must be solved.
            for i=1:NodeNum
                for j=1:NodeNum
                    
                    coeffs=zeros(1,inverse_order+1);
                    coeffs(1) = -covar_mat(i,j); % first include the constant
                    for kk=1:inverse_order
                       coeffs(kk+1)= PCweights(i,kk+1)*PCweights(j,kk+1)*factorial(kk);
                    end
                    coeffs = fliplr(coeffs);
                    r = roots(coeffs);   % roots only accepts decending-order coeffs
                    
                    count=0;
                    for kk=1:length(r)
                        if isreal(r(kk)) && (abs(r(kk))-1)<0.01  % 0.01 is the tolerance
                           count=count+1;
                           sol(count)=r(kk); 
                        end
                    end
                    
                    if count==0
                        % count=0 means no reasonable root found
                        error('No reasonable root exist!!!'); 
                    elseif count==1
                        covar_gaussian(i,j)=sol;
                    else                    
                       % If this error comes up, you can comment this error message and select one
                       % the reasonable root to continue
                       error('two or more reasonable roots exist!!!'); 
                       %covar_gaussian(i,j)=sol(1);                      
                    end

                end
            end

            
    end


end