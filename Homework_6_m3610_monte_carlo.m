
r = 500; % initial supporters in Group I
n = 10000; % people in Group I
s = 250; % initial supporters in Group II

p = .5; % as in problem statement
q = .1; % as in problem statement

k = floor((log((s*(1-p)/(n-r))/(1-q)^2)))/(log((1-p)/(1-q))) - 1 % optimal number of mentions

tt = 1*(1:(1*k));      %all time instances when we change the population size

max_final = 0;
min_final = Inf;

trials = 3000;
m0 = r + s;        %init population size

mm_end = [];    %final population size observed in all trial 

averaged_mm = zeros(size(tt));

for j=1:trials
    
    m = m0;
    mm = [m];
    r_1 = r;
    s_1 = s;
    
    for t=tt(2:end);
        r_1 = r_1 + sum(rand(1, n - r_1) < p);
        s_1 = s_1 - sum(rand(1, s_1) < q);
        m = r_1 + s_1;  %updated/current population size
        mm = [mm, m];  %population change history 
    end
    
    averaged_mm = averaged_mm + mm;
    
    if (max_final < mm(end))
        max_final = mm(end);
    else
        if (min_final > mm(end))
            min_final = mm(end);
        end
    end
    
    mm_end = [mm_end, mm(end)];
        
end

figure;  histogram(mm_end)