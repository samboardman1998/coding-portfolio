P = zeros(1,1);
D = zeros(1,1);
A_bar = zeros(1,1);
A_bound = zeros(1,1);
for p_inc = 1:99
    for d_inc = 1:99
        p = .5 + .5*(p_inc/100);
        d = .5 + .5*(d_inc/100);
        if p > 1/(2*d)
            V = zeros(1,1001);
            for V_inc = 1:1001
                mu_a = 0 + 1*((V_inc-1)/1000);
                V(V_inc) = max(mu_a,1-mu_a);
            end
            update = 1;
            while update > .00000000001
                update = 0;
                for V_inc = 1:1001
                    mu_a = 0 + 1*((V_inc-1)/1000);
                    mu_a_alpha = (p*mu_a)/(p*mu_a+(1-p)*(1-mu_a));
                    mu_a_alpha_inc = round(mu_a_alpha*1000+1);
                    mu_a_beta = ((1-p)*mu_a)/(p*(1-mu_a)+(1-p)*mu_a);
                    mu_a_beta_inc = round(mu_a_beta*1000+1);
                    old_V = V(V_inc);
                    V(V_inc) = max([mu_a 1-mu_a d*(mu_a*(p*V(mu_a_alpha_inc) + (1-p)*V(mu_a_beta_inc)) + (1-mu_a)*(p*V(mu_a_beta_inc) + (1-p)*V(mu_a_alpha_inc)))]);
                    update = max(update,V(V_inc) - old_V);
                end
            end
            disp('new paramter values')
            p
            d
            a_inc = 1001;
            a_bar = .5;
            while a_inc > 50
                a = 0 + 1*((a_inc-1)/1000);
                if V(a_inc) - a > .00005
                    a_bar = a + 1/1000;
                    a_inc = 50;
                else
                    a_inc = a_inc - 1;
                end    
            end
            a_bar
            V
            P = [P,p];
            D = [D,d];
            A_bar = [A_bar,a_bar];
            a_bound = (d*p+p-1)/(2*p-1);
            A_bound = [A_bound,a_bound];
        end
    end
end
P(1) = [];
D(1) = [];
A_bar(1) = [];
A_bound(1) = [];
A_diff = A_bound - A_bar;