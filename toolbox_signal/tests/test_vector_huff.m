% test of hufmann coding by concatening several token

% probability of having 0
t = .12;
n = 4096*2;
x = (rand(n,1)>t)+1;

% entropy lower bound
p = [t 1-t];
e =  -sum(p.*log2(p));

% create a new vector by lifting
q = 3;
n1 = ceil(n/q)*q;
x1 = x;
x1(end+1:n1) = 1;
x1 = reshape(x1,[q n1/q]);
[Y,X] = meshgrid(1:n1/q,0:q-1);
x1 = sum( (x1-1) .* (2.^X), 1 )' + 1;

% generate probability table
P = p(:); p = p(:);
for i=1:q-1
    Pold = P;
    P = [];
    for i=1:length(p)
        P = [P; Pold*p(i)];
    end
end

% compute the tree
T = compute_hufftree(P);
% do the coding
y = perform_huffcoding(x1,T,+1);
% average number of bits
e1 = length(y)/length(x);

disp(['Entropy=' num2str(e) ', Huffman(block size ' num2str(q) ')=' num2str(e1)]);