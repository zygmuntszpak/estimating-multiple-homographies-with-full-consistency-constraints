function u = unit(i, n)
% If i is a scalar, returns a unit vector.
% Otherwise returns a matrix of indicators u such that u(i(j),j) = 1 for all j.

u = zeros(n, length(i));
for j = 1:length(i)
  u(i(j),j) = 1;
end
