function h = vech(m)
% h = vech(m)
% h is the column vector of elements on or below the main diagonal of m.
% m must be square.

r = repmat((1:rows(m))', 1, cols(m));
c = repmat(1:cols(m), rows(m), 1);
% c <= r is same as tril(ones(size(m)))
h = m(find(c <= r));
