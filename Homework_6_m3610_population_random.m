

r = 0.7;        %rate of reproduction (newborn per current animal per year)

n = 1;          %how often do we "compound"
p = r/n;        %rate of reproduction per "time step"

T = 20;         %number of years
deltaT = 1/n;   %time step

tt = deltaT*(0:(n*T));      %all time instances when we change the population size

close all; 
figure; hold on;

max_final = 0;
min_final = Inf;

trials = 40;
m0 = 40;        %init population size

mm_end = [];    %final population size observed in all trial 



averaged_mm = zeros(size(tt));

tic
for j=1:trials
    m = m0;
    mm = [m];
   
    for t=tt(2:end);
        m = m + round(normrnd(m*p, sqrt(m*p*(1-p))));  %updated/current population size
        mm = [mm, m];  %population change history
    end
    
   % h = plot(tt, mm, 'g');
    
end
toc


tic
for j=1:trials
    m = m0;
    mm = [m];
    for t=tt(2:end);
        m = m + sum(rand(1,m) < p);  %updated/current population size
        mm = [mm, m];  %population change history
    end

   % h = plot(tt, mm, 'b');
    
end
toc

plot_handles = [h];

deterministic_mm = mm(1)*(1+p).^(0:(n*T));
h = plot(tt, deterministic_mm, 'r', 'LineWidth', 2);
plot_handles = [plot_handles, h];




%make it look nice  
xlabel('time');
ylabel('population');
str  = strcat('Exponential Population Growth:  m_0=', num2str(m0), ',  T=', num2str(T));  
title(str);
% legend(plot_handles, '    Normal approximation', '    Original stochastic', '    deterministic');

% Uncomment the following if you want to see a histogram of final
% population sizes (over all trials).
%
% figure;  hist(mm_end, 20)