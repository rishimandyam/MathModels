%% Euler's method

%a.
df2=@(x,y) y*sin(11*pi*x/6);
h=.01;
x0=0;
y0=1;
xrange=[0 1];
eulersmethod(df2,x0,y0,xrange,h)
eulersalgo(df2,x0,y0,xrange,.1)
legend("euler's method","displaying algorithm")
title("Euler's method for part a")
xlabel("x")
ylabel("y")
ylim([1 1.5])

%% b
x0=0;
y0=[1 2];
h=.01;
xrange=[0 5];
figure
eulersmethod(@system,x0,y0,xrange,h)
legend("y(x)","z(x)")
xlabel("x")
ylabel("y")
ylim([-12 5])
title("Euler's method for part b")

%% c
figure
df4=@(x,y) 3.5*y + 3*x - exp(y);
x0=0;
y0=0;
xrange=[0 10];
t=zeros(1,10);
for i=1:10
    % tic toc gives the computational time
    tic
    h=i/100;
    eulersmethod(df4,x0,y0,xrange,h)
    t(i)=toc;
end
title("Euler's method for part c with different step sizes.")
xlabel("x")
ylabel("y")
times=array2table(t,'VariableNames',{'.01','.02','.03','.04','.05','.06','.07','.08','.09','.1'});
%% Runge Kutta 2nd order part a
figure
df2=@(x,y) y*sin(11*pi*x/6);
h=.01;
x0=0;
y0=1;
xrange=[0 1];
rungekutta2(df2,x0,y0,xrange,h)
rk2algo(df2,x0,y0,xrange,.1)
title("Second order Runge Kutta for part a")
legend("Runge Kutta second order","Displaying algorithim")
xlabel("x")
ylabel("y")
ylim([1 1.5])

%% Runge Kutta 2nd order part b
figure
x0=0;
y0=[1 2];
h=.01;
xrange=[0 5];
rungekutta2(@system,x0,y0,xrange,h)
legend("y(x)","z(x)")
title("Second order Runge Kutta for part b")
ylim([-12 5])
xlabel("x")
ylabel("y")

%% Runge Kutta 2nd order part c
figure
df4=@(x,y) 3.5*y + 3*x - exp(y);
x0=0;
y0=0;
xrange=[0 10];
t=zeros(1,10);
for i=1:10
    tic
    h=i/100;
    rungekutta2(df4,x0,y0,xrange,h)
    t(i)=toc;
end
title("Runge Kutta 2nd order for part c with different step sizes.")
xlabel("x")
ylabel("y")
times=array2table(t,'VariableNames',{'.01','.02','.03','.04','.05','.06','.07','.08','.09','.1'});

%% Runge Kutta 4th order part a
figure
df2=@(x,y) y*sin(11*pi*x/6);
h=.01;
x0=0;
y0=1;
xrange=[0 1];
rungekutta4(df2,x0,y0,xrange,h)
rk4algo(df2,x0,y0,xrange,.1)
title("Fourth order Runge Kutta for part a")
legend("y(x)","Displaying algorithm")
xlabel("x")
ylabel("y")
ylim([1 1.5])

%%  Runge Kutta 4th order part b
figure
x0=0;
y0=[1 2];
h=.01;
xrange=[0 5];
rungekutta4(@system,x0,y0,xrange,h)
legend("y(x)","z(x)")
title("Fourth order Runge Kutta for part b")
ylim([-12 5])
xlabel("x")
ylabel("y")

%% Runge Kutta 4th order part c
figure
df4=@(x,y) 3.5*y + 3*x - exp(y);
x0=0;
y0=0;
xrange=[0 10];
t=zeros(1,10);
for i=1:10
    tic
    h=i/100;
    rungekutta4(df4,x0,y0,xrange,h)
    t(i)=toc;
end
title("Runge Kutta 4th order for part c with different step sizes.")
xlabel("x")
ylabel("y")
times=array2table(t,'VariableNames',{'.01','.02','.03','.04','.05','.06','.07','.08','.09','.1'});

%% subfn
function f=system(x,y_var)
% This subfunction defines the system of differential equations for Euler's
% method part b.
y=y_var(1);
z=y_var(2);
f(1)=-z;
f(2)=sin(3*x) + cos(2*x);
end

%% Euler's method
function eulersmethod(df,x0,y0,xrange,h)
% This function performs Euler's method with a given step size to calculate points that estimate
% the solution curve of a differential equation. This function graphs these
% points as a smooth curve over a given domain. This function is also
% capable of solving a system of differential equations and plotting the
% solutions.
% This function works by starting with an initial condition and using the 
% slope at the initial point to linearly approximate the next point by 
% using a linear model in the form of y=mx+b.

% This function requires five inputs.
% df is the differential equation which should be defined as an anonymous
% function or in the case of a system of differential equations, should be
% defined as a subfunction where the output is a vector that contains the
% differential equations. See the subfunction "system" to see how to define a system of differential equations in a subfunction. 
% The differential equation must have two and only two inputs, the second of 
% which can be a vector in cases where there are more than two variables.
% x0 is the independent variabel value for initial condition of the
% solution to the differential equation.
% y0 is the dependent variable value for initial condition of the solution to the differential
% equation. It is the value of the solution at x0. It
% should be input either as scalar for one differential equation or a
% vector when there is a system of equations.
% Example for y0:
% when x0 = 0
% a system of two differential equations with y(0)=1 and z(0)=3:
% y0 = [1 3]
% xrange is the x range which the solution(s) will be plotted over. It
% should be a two element row vector where the first element is the start
% of the range and the second element is the end of the range. The start of
% the range can not be less than the x0 value. For example, if x0 =0, the
% x range can not start at -1, it must start at 0 or greater.
% Example for xrange:
% where x0 = 0
% plotting from x=0 to x=10:
% xrange = [0 10]
% h is the step size for Euler's method. This means that this is the amount
% the dependent variable value changes for each new independent variable value approximation that is calculated. 
% h must a positive scalar value. The closer h is to zero, the better that
% the approximate solution will be. Generally, an h of around .01
% works well.

steps=(xrange(2)-x0)/h;
size=length(y0);

% This if statement handles if the x range is not evenly divisible by the
% step size.
if rem(steps, 1)==0
    steps=steps;
    x=x0:h:xrange(2);
elseif rem(steps, 1)~=0
    steps=ceil(steps);
    x=[x0:h:xrange(2) xrange(2)];
end
y=zeros(size,steps);
for j=1:size
for i=1:steps
    y(:,1)=y0;
    % This calculates the next y value.
    y(:,i+1)=y(:,i) + (x(i+1)-x(i))*df(x(i),y(:,i))';
    % This accounts for if the solution starts approaching negative or
    % positive infinity.
    if abs(y(:,i+1)) > 100000
        disp('Error, the function appears to diverge.')
        return
    end
end
% These four lines of code adjust for if the x range does not start at the
% initial condition.
newx=find(x>=xrange(1));
newx=x(:,newx(1):end);
extra=length(x)-length(newx);
newy=y(:,extra+1:end);

% This for loop assigns different colors to each line.
for i=1:size
    if i==1
        color='k';
    elseif i==2
        color='b';
    elseif i==3
        color='r';
    end
    plot(newx,newy(i,:),'-','Color',color)
    hold on
end
end
end

%% Euler's method displaying algorithm
function eulersalgo(df,x0,y0,xrange,h)
% This works indentically to the Euler's method function except it
% plots a curve with markers after each step of Euler's method. This
% function is meant to be used with a larger step size in order to display the steps of Euler's method more
% clearly. For example, if a step size of .01 gives a smooth curve, a step
% size of .1 should be used in the function so that the straight lines
% for each step can be clearly seen.
steps=(xrange(2)-x0)/h;

if rem(steps, 1)==0
    steps=steps;
    x=x0:h:xrange(2);
elseif rem(steps, 1)~=0
    steps=ceil(steps);
    x=[x0:h:xrange(2) xrange(2)];
end

y=zeros(1,steps);
for i=1:steps
    y(1)=y0;
    y(i+1)=y(i) + h*df(x(i),y(i));
    if abs(y(i+1)) > 100000
        disp('Error, the function appears to diverge.')
        return
    end
end
newx=find(x>=xrange(1));
newx=x(:,newx(1):end);
extra=length(x)-length(newx);
newy=y(:,extra+1:end);
plot(newx,newy,'r-')
hold on
plot(newx,newy,'ko')
end

%% Runge Kutta second order
function rungekutta2(df,x0,y0,xrange,h)
% This function performs the second order Runge Kutta method with a given step size to calculate points that estimate
% the solution curve of a differential equation. This function graphs these
% points as a smooth curve over a given domain. This function is also
% capable of solving a system of differential equations and plotting the
% solutions.
% This function works by starting with an initial condition and using the 
% slope at the initial point, k1, to linearly approximate the next hypothetical
% point. At this next hypothetical point the slope, k2, is calculated. The two
% slopes are averaged and from the initial point, the next actual point is
% approximated using a linear model in the form of y=mx+b.

% This function requires five inputs.
% df is the differential equation which should be defined as an anonymous
% function or in the case of a system of differential equations, should be
% defined as a subfunction where the output is a vector that contains the
% differential equations. See the subfunction "system" to see how to define a system of differential equations in a subfunction. 
% The differential equation must have two and only two inputs, the second of 
% which can be a vector in cases where there are more than two variables.
% x0 is the independent variabel value for initial condition of the
% solution to the differential equation.
% y0 is the dependent variable value for initial condition of the solution to the differential
% equation. It is the value of the solution at x0. It
% should be input either as scalar for one differential equation or a
% vector when there is a system of equations.
% Example for y0:
% when x0 = 0
% a system of two differential equations with y(0)=1 and z(0)=3:
% y0 = [1 3]
% xrange is the x range which the solution(s) will be plotted over. It
% should be a two element row vector where the first element is the start
% of the range and the second element is the end of the range. The start of
% the range can not be less than the x0 value. For example, if x0 =0, the
% x range can not start at -1, it must start at 0 or greater.
% Example for xrange:
% where x0 = 0
% plotting from x=0 to x=10:
% xrange = [0 10]
% h is the step size for the second order Runge Kutta method. This means that this is the amount
% the dependent variable value changes for each new independent variable value approximation that is calculated. 
% h must a positive scalar value. The closer h is to zero, the better that
% the approximate solution will be. Generally, an h of around .01
% works well.

steps=(xrange(2)-x0)/h;
size=length(y0);

if rem(steps, 1)==0
    steps=steps;
    x=x0:h:xrange(2);
elseif rem(steps, 1)~=0
    steps=ceil(steps);
    x=[x0:h:xrange(2) xrange(2)];
end
y=zeros(size,steps);
k1=zeros(size,steps);
k2=zeros(size,steps);
for j=1:size
for i=1:steps
    y(:,1)=y0;
    k1(:,i)=df(x(i),y(:,i));
    k2(:,i)=df(x(i+1),y(:,i) + k1(:,i)*(x(i+1)-x(i)));
    y(:,i+1)=y(:,i) + (x(i+1)-x(i))*((k1(:,i)+k2(:,i))/2);
    if abs(y(:,i+1)) > 100000
        disp('Error, the function appears to diverge.')
        return
    end
end
newx=find(x>=xrange(1));
newx=x(:,newx(1):end);
extra=length(x)-length(newx);
newy=y(:,extra+1:end);
for i=1:size
    if i==1
        color='k';
    elseif i==2
        color='b';
    elseif i==3
        color='r';
    end
    plot(newx,newy(i,:),'-','Color',color)
    hold on
end
end
end

%% Runge Kutta second order displaying algorithm
function rk2algo(df,x0,y0,xrange,h)
% This works indentically to the second order Runge Kutta method function except it
% plots a curve with markers after each step of the second order Runge Kutta method. This
% function is meant to be used with a larger step size in order to display the steps of Euler's method more
% clearly. For example, if a step size of .01 gives a smooth curve, a step
% size of .1 should be used in the function so that the straight lines
% for each step can be clearly seen.
steps=(xrange(2)-x0)/h;
size=length(y0);

if rem(steps, 1)==0
    steps=steps;
    x=x0:h:xrange(2);
elseif rem(steps, 1)~=0
    steps=ceil(steps);
    x=[x0:h:xrange(2) xrange(2)];
end
y=zeros(size,steps);
k1=zeros(size,steps);
k2=zeros(size,steps);
for j=1:size
for i=1:steps
    y(:,1)=y0;
    k1(:,i)=df(x(i),y(:,i));
    k2(:,i)=df(x(i+1),y(:,i) + k1(:,i)*(x(i+1)-x(i)));
    y(:,i+1)=y(:,i) + (x(i+1)-x(i))*((k1(:,i)+k2(:,i))/2);
    if abs(y(:,i+1)) > 100000
        disp('Error, the function appears to diverge.')
        return
    end
end
newx=find(x>=xrange(1));
newx=x(:,newx(1):end);
extra=length(x)-length(newx);
newy=y(:,extra+1:end);
plot(newx,newy,'r-')
hold on
plot(newx,newy,'ko')
end
end
%% Runge Kutta fourth order
function rungekutta4(df,x0,y0,xrange,h)
% This function performs the fourth order Runge Kutta method with a given step size to calculate points that estimate
% the solution curve of a differential equation. This function graphs these
% points as a smooth curve over a given domain. This function is also
% capable of solving a system of differential equations and plotting the
% solutions.
% This function works by starting with an initial condition and using the 
% slope at the initial point, k1, to linearly approximate the next hypothetical
% point which is half way to the next actual point. At this next hypothetical point the slope, k2, is calculated. 
% Then by using k2, the hypothetical point from the first step is
% recalculated and at this point the slope, k3, is calculated. Then, by
% using k3, a new hypothetical point is calculated that is the same x
% distance away from the initial point as the next actual point. At this
% new hypothetical point, the slope, k4, is calculated
% The four slopes are combined into a weighted averaged with k2 and k3
% having twice the weight of k1 and k4.
% From the initial point, the next actual point is approximated using a linear model in the form of y=mx+b.

% This function requires five inputs.
% df is the differential equation which should be defined as an anonymous
% function or in the case of a system of differential equations, should be
% defined as a subfunction where the output is a vector that contains the
% differential equations. See the subfunction "system" to see how to define a system of differential equations in a subfunction. 
% The differential equation must have two and only two inputs, the second of 
% which can be a vector in cases where there are more than two variables.
% x0 is the independent variabel value for initial condition of the
% solution to the differential equation.
% y0 is the dependent variable value for initial condition of the solution to the differential
% equation. It is the value of the solution at x0. It
% should be input either as scalar for one differential equation or a
% vector when there is a system of equations.
% Example for y0:
% when x0 = 0
% a system of two differential equations with y(0)=1 and z(0)=3:
% y0 = [1 3]
% xrange is the x range which the solution(s) will be plotted over. It
% should be a two element row vector where the first element is the start
% of the range and the second element is the end of the range. The start of
% the range can not be less than the x0 value. For example, if x0 =0, the
% x range can not start at -1, it must start at 0 or greater.
% Example for xrange:
% where x0 = 0
% plotting from x=0 to x=10:
% xrange = [0 10]
% h is the step size for the fourth order Runge Kutta method. This means that this is the amount
% the dependent variable value changes for each new independent variable value approximation that is calculated. 
% h must a positive scalar value. The closer h is to zero, the better that
% the approximate solution will be. Generally, an h of around .01
% works well.
steps=(xrange(2)-x0)/h;
size=length(y0);

if rem(steps, 1)==0
    steps=steps;
    x=x0:h:xrange(2);
elseif rem(steps, 1)~=0
    steps=ceil(steps);
    x=[x0:h:xrange(2) xrange(2)];
end
y=zeros(size,steps);
k1=zeros(size,steps);
k2=zeros(size,steps);
k3=zeros(size,steps);
k4=zeros(size,steps);

for j=1:size
for i=1:steps
    y(:,1)=y0;
    k1(:,i)=df(x(i),y(:,i));
    k2(:,i)=df(((x(i+1)-x(i))/2)+x(i),y(:,i) + k1(:,i)*((x(i+1)-x(i))/2));
    k3(:,i)=df(((x(i+1)-x(i))/2)+x(i),y(:,i) + k2(:,i)*((x(i+1)-x(i))/2));
    k4(:,i)=df(x(i+1),y(:,i) + k3(:,i)*(x(i+1)-x(i)));
    y(:,i+1)=y(:,i) + (x(i+1)-x(i))*((k1(:,i) + 2.*k2(:,i) + 2.*k3(:,i) + k4(:,i))/6);
    if abs(y(:,i+1)) > 100000
        disp('Error, the function appears to diverge.')
        return
    end
end
newx=find(x>=xrange(1));
newx=x(:,newx(1):end);
extra=length(x)-length(newx);
newy=y(:,extra+1:end);
for i=1:size
    if i==1
        color='k';
    elseif i==2
        color='b';
    elseif i==3
        color='r';
    end
    plot(newx,newy(i,:),'-','Color',color)
    hold on
end
end
end
%% Runge Kutta fourth order displaying algorithm
function rk4algo(df,x0,y0,xrange,h)
% This works indentically to the fourth order Runge Kutta method function except it
% plots a curve with markers after each step of the fourth order Runge Kutta method. This
% function is meant to be used with a larger step size in order to display the steps of Euler's method more
% clearly. For example, if a step size of .01 gives a smooth curve, a step
% size of .1 should be used in the function so that the straight lines
% for each step can be clearly seen.
steps=(xrange(2)-x0)/h;
size=length(y0);

if rem(steps, 1)==0
    steps=steps;
    x=x0:h:xrange(2);
elseif rem(steps, 1)~=0
    steps=ceil(steps);
    x=[x0:h:xrange(2) xrange(2)];
end
y=zeros(size,steps);
k1=zeros(size,steps);
k2=zeros(size,steps);
k3=zeros(size,steps);
k4=zeros(size,steps);

for j=1:size
for i=1:steps
    y(:,1)=y0;
    k1(:,i)=df(x(i),y(:,i));
    k2(:,i)=df(((x(i+1)-x(i))/2)+x(i),y(:,i) + k1(:,i)*((x(i+1)-x(i))/2));
    k3(:,i)=df(((x(i+1)-x(i))/2)+x(i),y(:,i) + k2(:,i)*((x(i+1)-x(i))/2));
    k4(:,i)=df(x(i+1),y(:,i) + k3(:,i)*(x(i+1)-x(i)));
    y(:,i+1)=y(:,i) + (x(i+1)-x(i))*((k1(:,i) + 2.*k2(:,i) + 2.*k3(:,i) + k4(:,i))/6);
    if abs(y(:,i+1)) > 100000
        disp('Error, the function appears to diverge.')
        return
    end
end
newx=find(x>=xrange(1));
newx=x(:,newx(1):end);
extra=length(x)-length(newx);
newy=y(:,extra+1:end);
plot(newx,newy,'r-')
hold on
plot(newx,newy,'ko')
end
end
