F = @(t,y) 0 ; ode23(F,[0 10],1)

F = @(t,y) log(y) ; ode23(F,[0 10],1) 
F = @(t,y) 1/y ; ode23(F,[0 10],7,1.e-6)

F = @(t,y) [y(2); -y(1)];
ode23(F,[0 2*pi],[1; 0])

opts = odeset('reltol',1.e-4,'abstol',1.e-6, ...
'outputfcn',@odephas2);
ode23(F,[0 2*pi],[1; 0],opts)

opts = odeset('reltol',1.e-4,'abstol',5.e-6, ...
'outputfcn',@odephas2);
ode23(F,[0 2*pi],[1; 0],opts)

opts = odeset('reltol',1.e-4,'abstol',1.e-9, ...
'outputfcn',@odephas2);
ode23(F,[0 2*pi],[1; 0],opts)

ode23(@twobody,[0 2*pi],[1; 0; 0; 1])

y0 = [1; 0; 0; .25];
ode23(@twobody,[0 2*pi],y0)

y0 = [1; 0; 0; .25];
[t,y] = ode23(@twobody,[0 10*pi],y0,1.e-20)
plot(y(:,1),y(:,2),'-',0,0,'ro')
axis equal

