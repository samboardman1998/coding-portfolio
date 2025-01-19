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

mean_in_degree = E/V % 3.6773
mean_out_degree = E/V % 3.6773
in_degree = sum(A);
out_degree = sum(A,2);
% histogram(in_degree); % exponential
% histogram(out_degree); % {1,...,10}-valued; median: 1; mode: 1; next frequent: 10; heavy-tailed distribution 

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

T = 500;

s = 0.85;
s_add = ((1-s)/V)*ones(V,V);
N_tilde = s*N + s_add;
% N_tilde = sparse(N_tilde);
N_tilde_T = transpose(N_tilde);
% N_tilde_T = sparse(N_tilde_T);
r_0 = ones(V,1);
r_0 = (1/V)*r_0;
r_0 = (N_tilde_T)^T*r_0;
r_0 = transpose(r_0);

r_1 = ones(1,V);
r_1 = (1/V)*r_1;
prod = r_1;
for n = 1:T-1
   prod = prod*N_tilde;
   r_1 = r_1 + prod; 
end
r_1 = (1/T)*r_1;

r_2 = zeros(1,V);
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

r_3 = zeros(1,V);
for i = 1:V
    for t = 1:T
        tau = geornd(s) + 1;
        X = zeros(1,tau);
        X_1 = i;
        X(1,1) = X_1;
        total = 0;
        count = 0;
        for k = 2:tau
            ran = rand;
            if ran < 1-s
                X(1,k) = i;
            else
                while total < ran
                    count = count + 1;
                    total = total + N_tilde(X(1,k-1),count);
                end
                X(1,k) = count;
                total = 0;
                count = 0;
            end
        end
        r_3(1,X(1,tau)) = r_3(1,X(1,tau)) + 1;
    end
end
r_3 = (1/(T*V))*r_3;

r_4 = zeros(1,V);
for t = 1:T
    tau = geornd(1-s) + 1;
    X = zeros(1,tau);
    X_1 = unidrnd(V);
    X(1,1) = X_1;
    total = 0;
    count = 0;
    r_4(1,X(1,1)) = r_4(1,X(1,1)) + 1;
    for k = 2:tau
        ran = rand;
        while total < ran
            count = count + 1;
            total = total + N_tilde(X(1,k-1),count);
        end
        X(1,k) = count;
        total = 0;
        count = 0;
        r_4(1,X(1,k)) = r_4(1,X(1,k)) + 1;
    end
end
r_4 = ((1-s)/T)*r_4;

