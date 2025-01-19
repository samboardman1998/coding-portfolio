Edges = readmatrix('p2p-Gnutella04.txt');
E_0 = size(Edges);
E = E_0(1,1);
V = max(max(Edges));

A = zeros(V,V);
for i = 1:E
    j_1 = Edges(i,1);
    j_2 = Edges(i,2);
    A(j_1,j_2) = 1;
end

mean_in_degree = E/V
mean_out_degree = E/V
in_degree = sum(A);
out_degree = sum(A,2);
histogram(out_degree) % {1,...,10}-valued; heavy-tailed distribution 

N = zeros(V,V);
for i = 1:V
    if out_degree(i,1) > 0   
        for j = 1:V
            N(i,j) = A(i,j)/out_degree(i,1);
        end
    else
        N(i,i) = 1;
    end
end

T = 12000; % number of iterations to perform for each algorithm below

s = 0.85; % scale factor as in PageRank 
s_add = ((1-s)/V)*ones(V,V);
N_tilde = s*N + s_add; % scaled stochastic matrix
N_tilde_T = transpose(N_tilde);
r_0 = ones(V,1); % vector of PageRank values via power iteration
r_0 = (1/V)*r_0;
prod = r_0;
for n = 1:T
   prod = (N_tilde_T)*prod; 
end
r_0 = prod;
r_0 = transpose(r_0);

in_neighbor_out_degree = zeros(1,V);
for i = 1:V
    count = 0;
    total = 0;
    if in_degree(1,i) ~= 0
        for j = 1:V
            if A(j,i) == 1
                count = count + 1;
                total = total + out_degree(j,1);
            end
        end
        in_neighbor_out_degree(1,i) = total/count;
    else
        in_neighbor_out_degree(1,i) = 0;
    end
end

out_neighbor_in_degree = zeros(V,1);
for i = 1:V
    count = 0;
    total = 0;
    if out_degree(i,1) ~= 0
        for j = 1:V
            if A(i,j) == 1
                count = count + 1;
                total = total + in_degree(1,j);
            end
        end
        out_neighbor_in_degree(i,1) = total/count;
    else
        out_neighbor_in_degree(i,1) = 0;
    end
end

x = in_neighbor_out_degree;
y = r_0;
scatter(x,y);

r_1 = ones(1,V); % vector of PageRank values via Monte Carlo Algorithm 1 (expected number of visits by Markov chain) 
r_1 = (1/V)*r_1;
prod = r_1;
for n = 1:T-1
   prod = prod*N_tilde;
   r_1 = r_1 + prod; 
end

r_1 = (1/T)*r_1;

r_2 = zeros(1,V); % vector of PageRank values via Monte Carlo Algorithm 2 (simulated number of visits by Markov chain)  
for t = 1:T
    tau = geornd(s) + 1;
    X = zeros(1,tau);
    X_1 = unidrnd(V);
    X(1,1) = X_1;
    total = 0;
    count = 0;
    for k = 2:tau
        ran = rand;
        while total < ran
            count = count + 1;
            total = total + N_tilde(X(1,k-1),count);
        end
        X(1,k) = count;
        total = 0;
        count = 0;
    end
    r_2(1,X(1,tau)) = r_2(1,X(1,tau)) + 1;
end
r_2 = (1/T)*r_2; 


tv_1 = 0; % total-variation distance of Monte Carlo Algorithm 1 PageRank vector from the PageRank vector via power iteration
tv_2 = 0; % total-variation distance of Monte Carlo Algorithm 2 PageRank vector from the PageRank vector via power iteration

for i = 1:V
    tv_1 = tv_1 + abs(r_1(1,i)-r_0(1,i));
    tv_2 = tv_2 + abs(r_2(1,i)-r_0(1,i));
end

tv_1 = tv_1/2 
tv_2 = tv_2/2