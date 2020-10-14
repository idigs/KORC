function valout = maxval(valin)
%VAI 2-21-06 produces largest single value in array
valout = max(reshape(valin,prod(size(valin)),1));