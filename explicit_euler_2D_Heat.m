% Domain and Grid initialization
domain = 1;
n = 41;
hx = domain/(n-1);
hy = domain/(n-1);
gridx= (1/hx)+1;
gridy= (1/hy)+1;
Hx= 0:hx:1;
Hy= 1:-hy:0; 

% Boundary conditions
w= zeros(gridx,gridy);
w(1,:)= 0;     
w(gridx,:)=1;  
for i= gridy:-1:1
    w(i,1) = 1-(Hy(i))^3;
    w(i,end) = 1-sin(pi*Hy(i)/2);
end

% Time Parameters
dt = 0.0001;     
total_t = 0.16;  
t_n = total_t/dt; 
a= dt/(2*hx*hx);
b= dt/(2*hy*hy);

% Error and Iteration Estimation
errStable = 1e-4;
errValue=100;
w_current= w;
count=1;

% Used forQuestion3
x_coord= find(Hy==0.4);
y_coord= find(Hx==0.4);

% w calculation 
for timestep = 1:t_n  
     for i= 2:gridx-1
       for j= 2:gridy-1
         w(i,j)= a*(w_current(i-1,j) + w_current(i+1,j) -2*w_current(i,j)) +  b*(w_current(i,j-1) + w_current(i,j+1) - 2*w_current(i,j)) + w_current(i,j);
       end    
     end

     % Error and Iteration Update
     errValue= max(max(abs(w_current-w)));
     w_current = w;
     count = count +1;
     tc= timestep*dt;
     
     % Required Plots
     % Question 4
     if(tc ==0.01 || tc ==0.02 || tc ==0.04 || tc ==0.08 || tc ==0.16 )
        figure
        contourf(Hx,Hy,w);
        colorbar;
        xlabel('x');
        ylabel('y')
        title('temp contours at time=',tc);
        hold off
     end

     % Question 3
      if(mod(timestep,10)==0)
        figure(6)
        plot(tc,w(y_coord,x_coord),'.')
        hold on;
        xlabel('time');
        ylabel('temp');
        axis([0 0.16 0 0.4]);
        title('temperature variation w.r.t time at x=y=0.4, dt = 0.0001')
      end
 
      if(tc==0.16)
       hold off;
       figure(7)
       plot(w(:,y_coord),Hy,'-')
       xlabel('y ');
       ylabel('temp variation along y-axis at x=0.4 and t=0.16');
       title('temp variation at x=0.4, dt = 0.0001 and t=0.16')
      end       
end

  
