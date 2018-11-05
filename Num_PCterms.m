function [P]=Num_PCterms(Dim,order)
%to calculate the number of PC terms for specified Dimension and order per
%eqn 3.54 in book by Roger Ghanem
P=1;
for i=1:order
    temp=1; 
    for j=0:(i-1)   
        temp=temp*(Dim+j);    
    end
    P=P+temp/factorial(i);           %number of PC terms with required dim and order
end

