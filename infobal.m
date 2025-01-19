Edges = readmatrix("edges.csv");
Edges = Edges + ones(size(Edges,1),size(Edges,2));
E_0 = size(Edges);
E = E_0(1,1);
V = max(max(Edges));

A = zeros(V,V);
for i = 1:E
    j_1 = Edges(i,1);
    j_2 = Edges(i,2);
    A(j_1,j_2) = 1;
end

G = digraph(A);

A_1 = sparse(A);

comp_list = conncomp(G); % for each node, its component

node_list = conncomp(G,"OutputForm","cell"); % for each component, its nodes

node_list_mat = cell2mat(node_list);

C = size(node_list,2); % number of components

l_2_list = [];
l_2reduc_list = [];
l_2base_list = [];
opt_packets = [];

for t = 1:10

bias_values = zeros(1,V); % for each node, its bias value

for i = 1:V
    bias_values(i) = randi(201) - 101; % random integer in [-100,100]
end

mean_bias_values = zeros(1,C); % the mean bias value for each component

for i = 1:C
    nodes = cell2mat(node_list(i));
    total = 0;

    for j = 1:size(nodes,2)
        total = total + bias_values(nodes(j));
    end

    mean_bias_values(i) = round(total/size(nodes,2)); % the mean bias value for the component, as used for the l_2^2 solution

end

connectivity = zeros(C,C); % (i,j)th entry is indicator for edge from component i to j (i < j)

for i = 1:C-1
    I = cell2mat(node_list(i));
    for j = i+1:C
        J = cell2mat(node_list(j));
        for v_i = 1:size(node_list(i),2)
            for v_j = 1:size(node_list(j),2)
                if distances(G,I(v_i),J(v_j)) <= V
                    connectivity(i,j) = 1;
                end
            end
        end
    end
end

x = zeros(1,C); % gives solution of packets for each component
packets = zeros(1,C); % |x|

for j = 1:C
    x(j) = -mean_bias_values(j);
    for i = 1:j-1
        if connectivity(i,j) == 1
            x(j) = x(j) - x(i);
        end
    end
    packets(j) = abs(x(j));

end

l_2 = 0;

for i = 1:V
    l_2 = l_2 + (bias_values(i) - mean_bias_values(comp_list(i)))^2;
end

k = 8717; % budget of packets
            
x_reduc = zeros(1,C);

reduc_factor = min(k/sum(packets),1);

for j = 1:C
    x_reduc(j) = fix(reduc_factor*x(j));
end

component_bias_shift = zeros(1,C);
reduc_bias_values = zeros(1,V);

for j = 1:C
    component_bias_shift(j) = x_reduc(j); 
    for i = 1:j-1
        if connectivity(i,j) == 1
            component_bias_shift(j) = component_bias_shift(j) + x_reduc(i);
        end
    end
end

for i = 1:V
    reduc_bias_values(i) = component_bias_shift(comp_list(i)) + bias_values(i);
end

l_2reduc = 0;

for i = 1:V
    l_2reduc = l_2reduc + (reduc_bias_values(i))^2;
end

l_2base = 0;

for i = 1:V
    l_2base = l_2base + (bias_values(i))^2;
end

l_2_list(t) = l_2;
l_2reduc_list(t) = l_2reduc;
l_2base_list(t) = l_2base;
opt_packets(t) = sum(packets);

end




% x_reduc vs. x_trunc

% quantify number of packets needed in this setting and then consider truncation approximation (with multiple datasets and perhaps derived stochastic block model)    

% node_list = zeros(V);
% augment = [];

% for i = 1:V
    % augment = [node_list(comp_list(i)),i];
    % node_list(comp_list(i)) = augment;
% end