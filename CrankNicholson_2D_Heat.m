% Domain and Grid Initialization
domain = 1;
n = 41;
hx = domain/(n-1);
hy = domain/(n-1);
gridx= (1/hx)+1;
gridy= (1/hy)+1;
Hx= 0:hx:1;
Hy= 1:-hy:0; 

% Boundary Conditions
w= zeros(gridx,gridy);
w(1,:)= 0;     
w(gridx,:)=1;
for i= gridy:-1:1
    w(i,1) = 1-(Hy(i))^3;
    w(i,end) = 1-sin(pi*Hy(i)/2);
end
  
% Time Parameters
dt = 0.001;     
total_t = 0.16;  
t_n = total_t/dt; 
a= dt/(2*hx*hx);
b= dt/(2*hy*hy);
alpha=a;

% Used for Question 3
x_coord= find(Hy==0.4);
y_coord= find(Hx==0.4);

% Co-efficient Matrix T construction
n=gridx-2;
m=n*n;
T= zeros(m);
for i=1:m
    T(i,i)= (1+4*alpha);
end

for i=1:m-n
    for j=i+n
        T(i,j)=-alpha;
    end
end

for i=n+1:m
    for j=i-n
        T(i,j)=-alpha;
    end
end

for i=1:m-1
    for j=i+1
        T(i,j)=-alpha;
    end
end

for i=1:m-1
    for j=i+1
       if(mod(i,n)==0)
            T(i,j)=0;
       end
    end
end

for i=2:m
    for j=i-1
        T(i,j)=-alpha;
    end
end

for i=2:m
    for j=i-1
        if(mod(i,n)==1)
            T(i,j)=0;
        end
    end
end

% Error and Iteration Estimation
errStable = 1e-5;
errValue= 10000;
w_current = w;
count=1;

for timestep = 1:t_n
  m=n*n;
  LB_var=n;
  RB_var=n;
  BB_var=m;
  n3 = n;
  R= ones(m,1);
     
    for i=3:gridx-2
      for j=3:gridx-2
          R(j-1+n*(i-2))=  alpha*(w_current(i+1,j) + w_current(i-1,j)+ w_current(i,j+1)+ w_current(i,j-1)-4*w_current(i,j)) + w_current(i,j);
        
      end
    end
 
 for i=3:gridx-2
      j=2;
         R(LB_var+1)=  alpha*(w_current(i+1,j) + w_current(i-1,j)+ w_current(i,j+1)+ w_current(i,j-1)-4*w_current(i,j)) + w_current(i,j)+ alpha*w(i,j-1);
     
      LB_var=LB_var+n;
      
 end
                                                                                           
 for i=3:gridx-2
      j=n+1;
         R(RB_var+n)=  alpha*(w_current(i+1,j) + w_current(i-1,j)+ w_current(i,j+1)+ w_current(i,j-1)-4*w_current(i,j)) + w_current(i,j)+ alpha*w_current(i,j+1);
     
      RB_var=RB_var+n;
 end
  
 for j=3:gridx-2
     i=2;
     R(j-1)=alpha*(w_current(i+1,j) + w_current(i-1,j)+ w_current(i,j+1)+ w_current(i,j-1)-4*w_current(i,j)) + w_current(i,j) + alpha*w_current(i-1,j);
 end
 
 for j=3:gridx-2
     for i=n+1
        R(BB_var-n+2)=alpha*(w_current(i+1,j) + w_current(i-1,j)+ w_current(i,j+1)+ w_current(i,j-1)-4*w_current(i,j)) + w_current(i,j) + alpha*w_current(i+1,j);
     end
     BB_var=BB_var+1;
 end

 R(m)= alpha*(w_current(n+2,n+1) + w_current(n,n+1)+ w_current(n+1,n+2)+ w_current(n+1,n)-4*w_current(n+1,n+1)) + w_current(n+1,n+1) + alpha*(w_current(n+1,n+2)+w_current(n+2,n+1));
 R(m+1-n) = alpha*(w_current(n+2,2) + w_current(n,2)+ w_current(n+1,3)+ w_current(n+1,1)-4*w_current(n+1,2)) + w_current(n+1,2) + alpha*(w_current(n+1,1)+w_current(n+2,2));
 R(1)=alpha*(w_current(3,2) + w_current(1,2)+ w_current(2,3)+ w_current(2,1)-4*w_current(2,2)) + w_current(2,2) + alpha*(w_current(1,2)+w_current(2,1));
 R(n)=alpha*(w_current(3,n+1) + w_current(1,n+1)+ w_current(2,n+2)+ w_current(2,n)-4*w_current(2,n+1)) + w_current(2,n+1) + alpha*(w_current(1,n+1) + w_current(2,n+2));

 % w matrix construction
     W=T\R;
     w(2:end-1,2:end-1)= (reshape(W,[gridx-2,gridy-2]))';
     w_current=w;
     tc= dt*timestep;
    % Required Plots
     % Question 4
        if(tc ==0.01 || tc ==0.02 || tc ==0.04 || tc ==0.08 || tc ==0.16)
        figure
        contourf(Hx,Hy,w);
        xlabel('x');
        ylabel('y');
        colorbar;
        title("temp contour at time =",tc);
        hold off
        end
     
        % Question 3
        figure(6)
        plot(tc,w(y_coord,x_coord),'.')
        hold on;
        xlabel('time');
        ylabel('temp');
        axis([0 0.16 0 0.4]);
        title('temperature variation w.r.t time at x=y=0.4, dt = 0.001')
 
      if(tc==0.16)
       hold off;
       figure(7)
       plot(w(:,y_coord),Hy,'-')
       xlabel('y ');
       ylabel('temp variation along y-axis at x=0.4 and t=0.16');
       title('temp variation at x=0.4, dt = 0.001 and t=0.16')
      end
     
     
end