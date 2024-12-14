%The goal is to find the stable age distribution
%and the rate of population growth based on 
%the known rates of birth/death for various age groups.


%birth-rate by age
% source:   http://www.cdc.gov/nchs/data/statab/t991x07.pdf
% data taken from the year 1999.
%br(i) is the chance that a woman of age (i-1) will give birth this year
br = 1e-3*[
0
0
0
0
0
0
0
0
0
0
0.9
0.9
0.9
0.9
0.9
28.7
28.7
28.7
80.3
80.3
111
111
111
111
111
117.8
117.8
117.8
117.8
117.8
89.6
89.6
89.6
89.6
89.6
38.3
38.3
38.3
38.3
38.3
7.4
7.4
7.4
7.4
7.4
0.4
0.4
0.4
0.4
0.4
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0];

%death-rate by age
% source:  http://www.ssa.gov/oact/STATS/table4c6.html
% data taken for females from the year 2009.
%dr(i) is the chance that a woman of age (i-1) will die this year
dr = [
0.005728
0.000373
0.000241
0.000186
0.00015
0.000133
0.000121
0.000112
0.000104
0.000098
0.000094
0.000098
0.000114
0.000143
0.000183
0.000229
0.000274
0.000314
0.000347
0.000374
0.000402
0.000431
0.000458
0.000482
0.000504
0.000527
0.000551
0.000575
0.000602
0.00063
0.000662
0.000699
0.000739
0.00078
0.000827
0.000879
0.000943
0.00102
0.001114
0.001224
0.001345
0.001477
0.001624
0.001789
0.001968
0.002161
0.002364
0.002578
0.0028
0.003032
0.003289
0.003559
0.003819
0.004059
0.004296
0.004556
0.004862
0.005222
0.005646
0.006136
0.006696
0.007315
0.007976
0.008676
0.009435
0.010298
0.011281
0.01237
0.013572
0.014908
0.01644
0.018162
0.020019
0.022003
0.024173
0.026706
0.029603
0.032718
0.036034
0.039683
0.043899
0.048807
0.054374
0.060661
0.067751
0.075729
0.084673
0.094645
0.105694
0.117853
0.131146
0.145585
0.161175
0.17791
0.195774
0.213849
0.231865
0.249525
0.266514
0.282504
0.299455
0.317422
0.336467
0.356655
0.378055
0.400738
0.424782
0.450269
0.477285
0.505922
0.536278
0.568454
0.602561
0.638715
0.677038
0.71766
0.76072
0.806363
0.851378
0.893947
];

%br = br*0.1;

brresidual = .1*br;
prevY = -inf;

minR_1 = 1;

for k=10:-1:1
   
    br = k*brresidual; 



%now let's form a Leslie matrix
A = zeros(120, 120);

%births
A(1,:)  =  br';
%ageing & dying
A = A + diag(1-dr(1:end-1), -1);
%assume the same rate of dying beyond 119 years 
A(120,120) = 1-dr(120);


x = rand(120,1);
 
   for i=1:200
    x=A*x;     %next year's population break down
    %normalize to make x(i) just a proportion of poupulation of age (i-1)
    x = x/sum(x);  
    
    %show the resulting age distribution every five years
   % if (mod(i,5) == 0)
   %    bar(x); 
   %    title( strcat('i = ', num2str(i)) );
   %    pause(0.2);
   % end
    end

%recover the same stationary age distribution as
%the eigenvector corresponding to an eigenvalue of maximum modulus
%
%the fact that this works is related to the Perron-Frobenius Thm
%(a simplified version to be covered in the next lecture)
[v,d]=eig(A);
[c,j] = max(abs(diag(d)));
y=v(:,j);
if(sum(y) >= prevY)
    prevY = sum(y);
    minY = .1*k;
end
sum2 = 0;
for t = 15:64
    sum2 = sum2 + y(t + 1);
end
if(sum2/sum(y) >= .6)
    minR_1 = .1*k;
end
   
hold on;
%plot(y/sum(y),'r', 'Linewidth', 2, 'color',rand(1,3));

%legend('r = .1', 'r = .2','r = .3', 'r = .4', 'r = .5', 'r = .6', 'r = .7', 'r = .8', 'r = .9', 'r = 1.0') 

end

minY %answer to a.
minR_1 % answer to b.


br = 10*brresidual; 

%now let's form a Leslie matrix
A = zeros(120, 120);

%births
A(1,:)  =  br';
%ageing & dying
A = A + diag(1-dr(1:end-1), -1);
%assume the same rate of dying beyond 119 years 
A(120,120) = 1-dr(120);


x = rand(120,1);
 
   for i=1:200
    x=A*x;     %next year's population break down
    %normalize to make x(i) just a proportion of poupulation of age (i-1)
    x = x/sum(x);  
    
    %show the resulting age distribution every five years
   % if (mod(i,5) == 0)
   %    bar(x); 
   %    title( strcat('i = ', num2str(i)) );
   %    pause(0.2);
   % end
    end

%recover the same stationary age distribution as
%the eigenvector corresponding to an eigenvalue of maximum modulus
%
%the fact that this works is related to the Perron-Frobenius Thm
%(a simplified version to be covered in the next lecture)
[v,d]=eig(A);
[c,j] = max(abs(diag(d)));
y=v(:,j);

CMF = [0];
PMF = [];
sumage = y(1);

for t = 1:119 %age
    CMFnew = (sumage/sum(y)).^200; 
    %CMFnew
    CMF = [CMF, CMFnew];
    
    PMFnew = CMF(t+1) - CMF(t);
    PMF = [PMF, PMFnew];
    
    sumage = sumage + y(t+1);
end
    PMFnew = 1 - sum(PMF);
    PMF = [PMF, PMFnew];
     
    %figure; histogram(PMF)
    
    bar(0:119, PMF);
 