%% A
clc
clear

f = @(x) 2.5.*sin(.25.*x.^2) + x + .7; % function
fd = @(x) 2.5*.5.*x.*cos(.25*x.^2) + 1; % Derivative
ig = 0;
t = .00000000001;
n_it = 10000;
xspan = [-5,5];
yspan = [-5,25];

Newts(f,fd,ig,t,n_it,xspan,yspan)
legend('Function','yaxis','New x','Tangent Lines','location','eastoutside')

%% B 

clc
clear

f1 = @(x) sin(2.*x);  % functions
f2 = @(x) (x-pi).^3 + pi/4;

fd1 = @(x) 2.*cos(2.*x); % derivatives
fd2 = @(x) 3.*(x-pi).^2;

ig = .5;
t = .000001;
n_it = 5;
xspan = [-5,5];
yspan = [-5,5];

f = {f1,f2}; % function vector
fd = {fd1,fd2}; % derivative function vector

Newts(f,fd,ig,t,n_it,xspan,yspan)
legend('Function 1','Function 2','x solution','location','eastoutside')
%% C
clc
clear

f = @(x) x.^7 - 2.*x.^6 - 30.*x.^5 + 40.*x.^4 + 229.*x.^3 - 198.*x.^2 -360.*x;
fd = @(x) 7*x^6 - 12*x.^5 - 150*x.^4 + 160*x.^3 + 687*x.^2 - 396*x - 360;
ig = [-50,.5,30];
t = .00000000001;
n_it = 10000;
xspan = [-6,6];
yspan = [-5000,5000];

Newts(f,fd,ig,t,n_it,xspan,yspan)
legend('Function','yaxis','New x','tangent lines','location','eastoutside')

%% D
clc
clear

f = @(x) x.^2;
fd = @(x) 2.*x;
ig = 6;
t = [1,.001,.000001,.00000001];
xspan = [-10,10];
yspan = [-10,10];
n_it = 10000;

Newts(f,fd,ig,t,n_it,xspan,yspan)
legend('Function','New x','Tangent Lines','location','eastoutside')

%% Newtons Method Function 

% To use this function, the first input should be a function or character vector of functions
% second input: derivative function or character vector of derivative functions
% third input: initial guess vector
% fourth input: tolerance value
% fifth input: maximum number of iterations you want to run function for
% sixth input: vector of x boundaries of the plot
% seventh input: vector of y boundaries of the plot

% this function will use Newton's method to solve for the roots of the
% given function using the function and derivative inputs along with an
% initial guess, tolerance level and number of iterations.

% The function will also output a plot with the function and tangent lines
% that help visualize how the function calculates the roots

function Newts(f,fd,ig,t,n_it,xspan,yspan)
    % To use this function, the first input should be a function or character vector of functions
    % second input: derivative function or character vector of derivative functions
    % third input: initial guess vector
    % fourth input: tolerance value
    % fifth input: maximum number of iterations you want to run function for
    % sixth input: vector of x boundaries of the plot
    % seventh input: vector of y boundaries of the plot

    % this function will use Newton's method to solve for the roots of the
    % given function using the function and derivative inputs along with an
    % initial guess, tolerance level and number of iterations.

    % The function will also output a plot with the function and tangent lines
    % that help visualize how the function calculates the roots
    
    if length(f) > 1 % code within this 'if block' is used to execute Newtons's method for two functions 
        fn1 = cell2mat(f(1));
        fn2 = cell2mat(f(2));
        fd1 = cell2mat(fd(1));
        fd2 = cell2mat(fd(2));
        
        fplot(fn1)
        hold on
        fplot(fn2)
        xlim(xspan)
        ylim(yspan)
        x = ig;
        f = [fn1(x);fn2(x)];
        fnd = [];
        n=1;
        xval = [x];
        while abs(f(2,n)-f(1,n)) > t % Computation for multiple function inputs
            newfunc = [fn1(x);fn2(x)];
            f = [f newfunc];
            nfnd = [fd1(x); fd2(x)];
            fnd = [fnd nfnd]; 
            b = f(1,n+1) - fnd(1,n)*x;
            y=@(x) fnd(1,n)*x + b - fn2(x);
            newx = fsolve(y,x);
            x = newx;
            n = n + 1;
            xval = [xval, x];
            newy = [fn1(x);fn2(x)];
            yval = [f newy];
            tab = table(xval',(yval(1,2:end))',(yval(2,2:end))','VariableNames',{'x','f1(x)','f2(x)'})
        end
        plot(x,fn1(x),'go') % plots solution
        fprintf('The calculated value of x is %f \n',x)
        
    else % this block is used to compute Newton's for one function
    fplot(f)
    xlim(xspan)
    ylim(yspan)
    hold on
    yline(0)
    
    xn = linspace(-6,6,500);
    n = 1;
    for ii = 1:length(ig) % accounts for multiple initial guesses
        x = ig(ii);
        y = f(x);
        for jj = 1:length(t) % accounts for multiple tolerances
            t1 = t(jj);
            subplot(length(t),1,jj)
            if length(t) > 1
                fplot(f)
            end
            hold on
            while abs(y) > t1 & n <= n_it  % Newton's computation
                x = x - (y/fd(x)); 
                y = f(x);
                yn = fd(x)*(xn - x) + f(x);
                hold on 
                plot(x,y,'r*')
                plot(xn,yn,'k-')
                n = n + 1;
            end
            if n >= n_it % accounts for divergence in calculation
                disp("ERROR: The Newton's Method calculation diverged")
            end
            plot(x,f(x),'go')
            fprintf('The calculated value of x is %f \n',x)
        end
    end
    end
end
