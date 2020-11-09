
V=random('exp',2,6,2)

dim = size(C,2)

P = (1:size(C,1))';
E = [] ;
CE = [];
assert(isempty(P) || (size(P,2) == 1))
assert( isempty(E) || (size(E,2) == 2))


np = numel(P)
ne = size(E,1)
m = np + ne
n = size(V, 1)
c = size(C,1)
s=repmat(V,[1,1,c])
ss=repmat(C,[1,1,n])
sss=permute(repmat(C,[1,1,n]),[3,2,1])
sssss=sum(repmat(V,[1,1,c]) - permute(repmat(C,[1,1,n]),[3,2,1]))
D = permute(sum((repmat(V,[1,1,c]) - permute(repmat(C,[1,1,n]),[3,2,1])).^2,2),[1,3,2])

[minD,Cv] = min(D)

if(~all(size(unique(Cv)) == size(Cv)))
    warning('Multiple control vertices snapped to the same domain vertex');
end

bc = repmat(NaN,[n m])
bc(Cv(P),:) = eye(np,m);
indices = 1:n;

b = indices(any(~isnan(bc),2))
bc = bc(b,:)
bc(isnan(bc)) = 0

bc(any(bc,2),:)  = bc(any(bc,2),:) ./ repmat(sum(bc(any(bc,2),:),2),1,size(bc,2))
