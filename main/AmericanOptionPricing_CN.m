clear all

%when underlying asset price is 10, the solution should be the middle point
%of the generated u vector, i.e. x = 0, V = strikePrice*exp(0) = strikePrice
assetPrice = 10;
strikePrice = 10;
interestRate = 0.05;
timeToMaturity = 0.5;
volatility = 0.20;
%for put option
flag = 0;

%setting dx = dt
dx = 0.05;
Nminus = -200;
Nplus = 200;

M = 10;
dt = timeToMaturity/M;

%reference solution
a = zeros(6,1);
for i = 1:6
    [AssetPrice,OptionValue] = binprice(assetPrice,strikePrice,interestRate,timeToMaturity,2^(-i)*0.01,volatility,flag);
    a(i) = OptionValue(1,1);
end
disp('reference price')
disp(a)

disp('Pricing American Options using PSOR with Implicit-FD')

time = cputime;
u0 = PSOR_CN(2*dx,2*dt,M/2,Nplus/2,Nminus/2,volatility,interestRate);
p0 = median(u0);
t0 = cputime - time;

time = cputime;
u1 = PSOR_CN(dx,dt,M,Nplus,Nminus,volatility,interestRate);
p1 = median(u1);
t1 = cputime - time;

time = cputime;
u2 = PSOR_CN(dx/2,dt/2,2*M,2*Nplus,2*Nminus,volatility,interestRate);
p2 = median(u2);
t2 = cputime - time;

time = cputime;
u3 = PSOR_CN(dx/4,dt/4,4*M,4*Nplus,4*Nminus,volatility,interestRate);
p3 = median(u3);
t3 = cputime - time;

time = cputime;
u4 = PSOR_CN(dx/8,dt/8,8*M,8*Nplus,8*Nminus,volatility,interestRate);
p4 = median(u4);
t4 = cputime - time;

time = cputime;
u5 = PSOR_CN(dx/16,dt/16,16*M,16*Nplus,16*Nminus,volatility,interestRate);
p5 = median(u5);
t5 = cputime - time;

r0 = (p0-p1)/(p1-p2);
r1 = (p1-p2)/(p2-p3);
r2 = (p2-p3)/(p3-p4);
r3 = (p3-p4)/(p4-p5);

PSORt = [t0;t1;t2;t3;t4;t5]
PSORprice = strikePrice*[p0;p1;p2;p3;p4;p5]
PSORratio = [r0;r1;r2;r3]

disp('Pricing American Options using Penalty Method with implicit FD')

%the last two parameters are expected tolerance and penalty L
time = cputime;
u0 = Penalty_Method_CN(2*dx,2*dt,M/2,Nplus/2,Nminus/2,volatility,interestRate,10^-8,10^8);
p0 = median(u0);
t0 = cputime - time;

time = cputime;
u1 = Penalty_Method_CN(dx,dt,M,Nplus,Nminus,volatility,interestRate,10^-8,10^8);
p1 = median(u1);
t1 = cputime - time;

time = cputime;
u2 = Penalty_Method_CN(dx/2,dt/2,2*M,2*Nplus,2*Nminus,volatility,interestRate,10^-8,10^8);
p2 = median(u2);
t2 = cputime - time;

time = cputime;
u3 = Penalty_Method_CN(dx/4,dt/4,4*M,4*Nplus,4*Nminus,volatility,interestRate,10^-8,10^8);
p3 = median(u3);
t3 = cputime - time;

time = cputime;
u4 = Penalty_Method_CN(dx/8,dt/8,8*M,8*Nplus,8*Nminus,volatility,interestRate,10^-8,10^8);
p4 = median(u4);
t4 = cputime - time;

time = cputime;
u5 = Penalty_Method_CN(dx/16,dt/16,16*M,16*Nplus,16*Nminus,volatility,interestRate,10^-8,10^8);
p5 = median(u5);
t5 = cputime - time;

r0 = (p0-p1)/(p1-p2);
r1 = (p1-p2)/(p2-p3);
r2 = (p2-p3)/(p3-p4);
r3 = (p3-p4)/(p4-p5);

PENALTYt = [t0;t1;t2;t3;t4;t5]
PENALTYprice = strikePrice*[p0;p1;p2;p3;p4;p5]
PENALTYratio = [r0;r1;r2;r3]

disp('Pricing American Options using Policy Iteration with Implicit-FD')

time = cputime;
u0 = Policy_Iteration_CN(2*dx,2*dt,M/2,Nplus/2,Nminus/2,volatility,interestRate);
p0 = median(u0);
t0 = cputime - time;

time = cputime;
u1 = Policy_Iteration_CN(dx,dt,M,Nplus,Nminus,volatility,interestRate);
p1 = median(u1);
t1 = cputime - time;

time = cputime;
u2 = Policy_Iteration_CN(dx/2,dt/2,2*M,2*Nplus,2*Nminus,volatility,interestRate);
p2 = median(u2);
t2 = cputime - time;

time = cputime;
u3 = Policy_Iteration_CN(dx/4,dt/4,4*M,4*Nplus,4*Nminus,volatility,interestRate);
p3 = median(u3);
t3 = cputime - time;

time = cputime;
u4 = Policy_Iteration_CN(dx/8,dt/8,8*M,8*Nplus,8*Nminus,volatility,interestRate);
p4 = median(u4);
t4 = cputime - time;

time = cputime;
u5 = Policy_Iteration_CN(dx/16,dt/16,16*M,16*Nplus,16*Nminus,volatility,interestRate);
p5 = median(u5);
t5 = cputime - time;

r0 = (p0-p1)/(p1-p2);
r1 = (p1-p2)/(p2-p3);
r2 = (p2-p3)/(p3-p4);
r3 = (p3-p4)/(p4-p5);

POLICYt = [t0;t1;t2;t3;t4;t5]
POLICYprice = strikePrice*[p0;p1;p2;p3;p4;p5]
POLICYratio = [r0;r1;r2;r3]


plot(log(PSORt),log(abs(a(6)-PSORprice)),'--o')
hold on
plot(log(PENALTYt),log(abs(a(6)-PENALTYprice)),'--o')
hold on
plot(log(POLICYt),log(abs(a(6)-POLICYprice)),'--o')
hold off
title('American Option Pricing under Crank-Nicolson schema')
xlabel('log(runtime)')
ylabel('log(error)')
legend('PSOR','Penalty','Policy')


