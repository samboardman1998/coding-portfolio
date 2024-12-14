

r = 0.7;        %rate of reproduction (newborn per current animal per year)
K = 200;        %carrying capacity of population

n = 1;          %how often do we "compound"
p = r/n;        %rate of reproduction per "time step"

T = 20;         %number of years
deltaT = 1/n;   %time step

tt = deltaT*(0:(n*T));      %all time instances when we change the population size

close all; 
figure; hold on;

max_final = 0;
min_final = Inf;

trials = 200;
m0 = 30;        %init population size

mm_end = [];    %final population size observed in all trial 



averaged_mm = zeros(size(tt));

for j=1:trials
    m = m0;
    mm = [m];
    for t=tt(2:end);
        m = m + round(sum(rand(1,m) < p)*(1 - m/K));  %updated/current population size
        mm = [mm, m];  %population change history
    end

    h = plot(tt, mm, 'b');
    
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

plot_handles = [h];

averaged_mm = averaged_mm / trials;
h = plot(tt, averaged_mm, 'k', 'LineWidth', 5);
plot_handles = [plot_handles, h];

deterministic_mm = K./(((K - m0)./m0)*exp(1).^(-r*(0:(n*T))) + 1);
h = plot(tt, deterministic_mm, 'r', 'LineWidth', 2);
plot_handles = [plot_handles, h];

approximation_diff_end = abs(averaged_mm(end) - deterministic_mm(end))

spread = (max_final - min_final)
relative_spread = spread / deterministic_mm(end)

%find the range of final population sizes once we exclude the top 10% and bottom 10% of simulations 
mm_end_without_outliers = prctile(mm_end, [10, 90]);  
spread_without_outliers = (mm_end_without_outliers(2) - mm_end_without_outliers(1))
relative_spread_without_outliers = spread_without_outliers / deterministic_mm(end)


%make it look nice  
xlabel('time');
ylabel('population');
str  = strcat('Logistic Population Growth:  m_0=', num2str(m0), ',  T=', num2str(T));  
title(str);
legend(plot_handles, '    stochastic', '    averaged', '    deterministic');

% Uncomment the following if you want to see a histogram of final
% population sizes (over all trials).
%
% figure;  hist(mm_end, 20)






r = 1;          %number of distinct conversations per hour
K = 15;         %number of individuals in the group
 
T = 20;         %number of hours
deltaT = 1/n;   %time step
 
tt = deltaT*(0:(n*T));      %all time instances when we change the population size
 
close all; 
figure; hold on;
 
max_final = 0;
min_final = Inf;
 
trials = 200;
m0 = 1;        %initially one person knows (or creates) the rumor
 
mm_end = [];    %final number of people that know the rumor observed in all trial 
 
 
 
averaged_mm = zeros(size(tt));
 
for j=1:trials
    m = m0;
    mm = [m];
    for t=tt(2:end);
        m = m + round(r*2*(m/K)*(1 - m/K));  %updated/current number of people that know the rumor 
        mm = [mm, m];  %number-of-people-that-know-rumor change history
    end
 
    h = plot(tt, mm, 'b');
    
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
 
plot_handles = [h];
 
averaged_mm = averaged_mm / trials;
h = plot(tt, averaged_mm, 'k', 'LineWidth', 5);
plot_handles = [plot_handles, h];
 
 
deterministic_mm = K./(((K - m0)./m0)*exp(1).^(-r*(0:(n*T))) + 1);
h = plot(tt, deterministic_mm, 'r', 'LineWidth', 2);
plot_handles = [plot_handles, h];
 
approximation_diff_end = abs(averaged_mm(end) - deterministic_mm(end))
 
spread = (max_final - min_final)
relative_spread = spread / deterministic_mm(end)
 
%find the range of final population sizes once we exclude the top 10% and bottom 10% of simulations 
mm_end_without_outliers = prctile(mm_end, [10, 90]);  
spread_without_outliers = (mm_end_without_outliers(2) - mm_end_without_outliers(1))
relative_spread_without_outliers = spread_without_outliers / deterministic_mm(end)
 
 
%make it look nice  
xlabel('time');
ylabel('population');
str  = strcat('Exponential Population Growth:  m_0=', num2str(m0), ', T=', num2str(T));  
title(str);
legend(plot_handles, 'stochastic', 'averaged', 'deterministic');
 
% Uncomment the following if you want to see a histogram of final
% population sizes (over all trials).
%
% figure;  hist(mm_end, 20)
 




