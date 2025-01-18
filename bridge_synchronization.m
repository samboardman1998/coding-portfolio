n = 200;
p_0 = .5;
m = 1600;

v_T = [];

for t = 1:4000
    
 l_0 = sum(rand(1,n) < p_0);
 r_0 = n - l_0;

 p_new = r_0/n;

 v = [l_0];
    
 T = 1600;
 
 if(l_0 == n || l_0 == 0)
     T = 0;
 end
 
 for i = 1:m
   
   if(mod(i, 2) == 1)
      
     r_new = n - sum(rand(1,n) < p_new);
     v = [v, r_new];
     p_new = r_new/n;
     
     if((r_new == 0 || r_new == n) && T == 1600)
        
         T = i;
         
     end
       
   else
    
     l_new = sum(rand(1,n) < p_new);  
     v = [v, l_new];
     p_new = (n - l_new)/n;
     
     if((l_new == 0 || l_new == n) && T == 1600)
        
         T = i;
         
     end
       
   end
       
 end
 
 %v = v/n;
 %plot(1:length(v), v);
 
 v_T = [v_T, T];
    
end
 
 E_T = mean(v_T) % estimated expected value of T
 figure;  histogram(v_T)

