Edges = readmatrix("irvine_social_network.txt");
% Edges = Edges + ones(size(Edges,1),size(Edges,2)); (for Gnutella06)
E_0 = size(Edges);
E = E_0(1,1);
V = max(max(Edges));

A = zeros(V,V);
for i = 1:E
    j_1 = Edges(i,1);
    j_2 = Edges(i,2);
    A(j_1,j_2) = 1;
end

M = .2*A; % probability of edges

G = digraph(A);

A_1 = sparse(A);

bias_values = zeros(1,V); % for each node, its bias value

%for i = 1:V
%    x = normrnd(0,33);
%    bias_values(i) = x;
%end

for i = 1:V
    x = trnd(1);
    bias_values(i) = x;
end

l_2base = 0;

for i = 1:V
    l_2base = l_2base + (bias_values(i))^2;
end

% try solution in expectation (can reuse simulation part), even though lacks theoretical guarantee

comp_list = conncomp(G); % for each node, its component

node_list = conncomp(G,"OutputForm","cell"); % for each component, its nodes

C = size(node_list,2); % number of components

ordered_nodes = cell2mat(node_list);

prob_estimate = estimate(G,V,M,100);

x = zeros(1,V);
packets = zeros(1,V);

for v = 1:V
    x(ordered_nodes(v)) = -bias_values(ordered_nodes(v));
    for u = 1:v-1
        x(ordered_nodes(v)) = x(ordered_nodes(v)) - prob_estimate(ordered_nodes(u),ordered_nodes(v))*x(ordered_nodes(u));
    end
    packets(ordered_nodes(v)) = abs(x(ordered_nodes(v)));
end

% reduced_bias_values = seed(12,1,G,V,M,bias_values);

val = 0;
new_bias = bias_values;

for v = 1:V
    
    if x(ordered_nodes(v)) > 0
        val = 1;
        for t = 1:packets(ordered_nodes(v))
            new_bias = seed(ordered_nodes(v),val,G,V,M,new_bias);
        end
    end

    if x(ordered_nodes(v)) < 0
        val = -1;
        for t = 1:packets(ordered_nodes(v))
            new_bias = seed(ordered_nodes(v),val,G,V,M,new_bias);
        end
    end
    
end

l_2reduc = 0;

for i = 1:V
    l_2reduc = l_2reduc + (new_bias(i))^2;
end

function [prob_estimate] = estimate(G,V,M,T)

prob_estimate = zeros(V,V);

for t = 1:T

for u = 1:V

isExposed = zeros(V);

isExposed(u) = 1;
toTransmitFrom = [u];

while size(toTransmitFrom,2) > 0
current = toTransmitFrom(1);
toTransmitFrom = toTransmitFrom(2:size(toTransmitFrom,2));
neighbors = nearest(G,current,1);
num_neighbors = size(nearest(G,current,1),1);
for n = 1:num_neighbors
    if isExposed(neighbors(n)) == 0
        if rand < M(current,neighbors(n)) 
            isExposed(neighbors(n)) = 1;
            toTransmitFrom = [toTransmitFrom,neighbors(n)];
            prob_estimate(u,neighbors(n)) = prob_estimate(u,neighbors(n)) + 1;
        end
    end
end

end

end

end

prob_estimate = (1/T)*prob_estimate;

end

function [reduc_bias_values] = seed(u,packet_value,G,V,M,bias_values)

reduc_bias_values = bias_values;

isExposed = zeros(V);

isExposed(u) = 1;
toTransmitFrom = [u];
reduc_bias_values(u) = reduc_bias_values(u) + packet_value;

while size(toTransmitFrom,2) > 0
current = toTransmitFrom(1);
toTransmitFrom = toTransmitFrom(2:size(toTransmitFrom,2));
neighbors = nearest(G,current,1);
num_neighbors = size(nearest(G,current,1),1);
for n = 1:num_neighbors
    if isExposed(neighbors(n)) == 0
        if rand < M(current,neighbors(n)) 
            isExposed(neighbors(n)) = 1;
            toTransmitFrom = [toTransmitFrom,neighbors(n)];
            reduc_bias_values(neighbors(n)) = reduc_bias_values(neighbors(n)) + packet_value;
        end
    end
end

end

end