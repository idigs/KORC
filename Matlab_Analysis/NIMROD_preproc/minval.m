function valout = minval(valin)
%VAI 2-21-06 produces largest single value in array
valout = min(reshape(valin,prod(size(valin)),1));