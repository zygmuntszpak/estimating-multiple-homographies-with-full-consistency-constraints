function [ N ] = create_matrix_N( n )
%CREATE_MATRIX_N Summary of this function goes here
%   Detailed explanation goes here
 N = [n(1) n(4) n(5); 
      n(4) n(2) n(6);
      n(5) n(6) n(3)];

end

