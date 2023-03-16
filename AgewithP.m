clear all
%%initial parameter
initialtime = 0;
timestep = 0.1;
finaltime = 100;
agestep = 0.5;
finalage = 60;
initialfunction = @(a)exp(-a);
%%Íø¸ñ»®·Ö
timenode = initialtime:timestep:finaltime;
spacenode = agestep:agestep:finalage;
spaceinterval = length(spacenode);
initialvalue=initialfunction(spacenode');

%%Example 2
    R_0 = 2;
    alpha=0.1;
    beta_0 = @(a)R_0*exp(-alpha*a);
    mu_0 =@(a)1-alpha;
    A1=20;
    A2=30;
    mu_1= zeros(1,length(spacenode));
    for i = 1:spaceinterval
       a=i*agestep;
       if a>=A1&& a<=A2
          mu_1(i)=3;
       else
           mu_1(i)=0;
       end
   end
   G = @(p)exp(-p);
   Q = @(p)p./(1+p);

%%Example 3
%    alpha=0.1;
%    beta_0 = @(a)a*exp(-alpha*a);
%     mu_0 = @(a)(1+a*alpha)/(1+a);
%     mu_1= @(a)0.*a;
%     G = @(p)(1+p)^(-2);
%     Q = @(p)0.*p;
    T(1)=(1+agestep*mu_0 (spacenode(1)))^(-1);
    R(1)= agestep*  beta_0(spacenode(1))*T(1);
for i=1:length(spacenode)-1
    T(i+1)=T(i)*(1+agestep*mu_0 (spacenode(i+1)))^(-1);
  R(i+1)=R(i)+ agestep*  beta_0(spacenode(i+1))*T(i+1);
end
R(end);

A = zeros(spaceinterval);
for i = 1:spaceinterval-1 
    A(i,i) = -1/agestep;
    A(i + 1,i) = 1/agestep;
end
A(spaceinterval,spaceinterval)=-1/agestep;
%%Matrix B
B1 = zeros(1,spaceinterval);
M0=zeros(spaceinterval);
M1=zeros(spaceinterval);
for i=1:spaceinterval
   B1(1,i)=beta_0(spacenode(i));
   M0(i,i)=mu_0(spacenode(i));
   M1(i,i)=mu_1(i);
end
B2 = zeros(spaceinterval-1,spaceinterval);
B=[B1;B2];
I=eye(spaceinterval);
p= zeros( spaceinterval, length(timenode));
P = zeros( 1,length(timenode));
p(:,1) =  initialvalue;
P(1) =  agestep * sum(p(:,1));

for i=1:length(timenode)-1
p(:,i+1) = (I- timestep*A+timestep*M0+timestep*Q(P(i))*M1)\ (I+timestep*G(P(i))* B)*p(:,i);
P(i+1) =  agestep * sum(p(:,i+1));  
end
figure (1)
plot(timenode,P,'LineWidth',2)
xlabel('t_n')
ylabel('P^n')
hold on
figure (2)
mesh(timenode,spacenode,p)
xlabel('t_n')
ylabel('a')
zlabel('P^n')
figure(3)
semilogy (spacenode,p(:,end-1))

    
    