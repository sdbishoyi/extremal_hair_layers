function [row_idx,col_idx]=find_value_arr(arr,value)

value_arr=zeros(size(arr))+value;

tol=min(min(abs(value_arr-arr)));

[row_idx,col_idx]=find( abs((arr-value))<=tol );

end