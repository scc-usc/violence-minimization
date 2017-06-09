function [ ninfec, immset ] = gen_graphs(pmat, x0, T, k, p, seeds, pid)
pmat = sparse(pmat);
nn =length(pmat);
pid = unique(pid);
lmat = zero_diag(pmat>0);
N = length(pid);
t1 = 1./(N - sum(lmat, 1));
%outmat = repmat(t1, [nn 1])';
outmat = sparse(length(pmat));
pidlog = sparse(length(pmat), 1);
pidlog(pid) = 1;

for j=1:length(pid)
   discon = lmat(:, pid(j))<1 & pidlog;
   outmat(discon, pid(j)) = t1(j);
   outmat(pid(j), pid(j)) = 0;
end
npmat = p*pmat' + (1-p)*outmat';
npmat_T = (1/T)*npmat;
M = npmat;
for i=2:T
    M = M*npmat;
    npmat_T = npmat_T + M/T;
end

x00 = x0;
x0(:) = 0;
x0(pid) = x00(pid);

if(isempty(seeds))
    coeffs = ones(1, nn)*npmat_T;
    coeffs(x0==0) = 0;
    [xx, idx] = sort(coeffs, 'descend');
    immset = idx(1:k);
else
    immset = seeds;
end

ninfec = zeros(k, 1);
for i=1:k
    xt = x0;
    xt(immset(1:i)) = 0;
    ninfec(i) = sum(npmat_T*xt);
end
sum(npmat_T*x0)
end

