%% A tridiagonal matrix with random elements
m = 5;
n = m;

Amd = rand(1,1,m);
Asub = rand(1,1,m-1);
Asup = rand(1,1,m-1);
A = blktridiag(Amd,Asub,Asup);

X = linsolve(full(A),ones(m,1));

mmwrite('sample.mtx',A,'Sample Tridiagonal matrix');

