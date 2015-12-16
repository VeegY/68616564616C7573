function [rowptr, colind, vals] = csrwrite(filename,A)
% Bringe die Matrix A in ein CSR-Format,
% das die Löserklassen parallel einlesen können.
% rowptr und colind sind 0-basiert

[indj,indi,vals] = find(A');
nentr = nnz(A);
n = length(A);

rowptr = zeros(n,1);
for i=1:n
rowptr(i) = find(indi==i,1)-1;
end

colind = indj-1;

dlmwrite(strcat(filename,'_vals.csr'),vals);
dlmwrite(strcat(filename,'_rowptr.csr'),rowptr);
dlmwrite(strcat(filename,'_colind.csr'),colind);

end
