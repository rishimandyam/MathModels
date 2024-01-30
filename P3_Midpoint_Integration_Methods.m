%% Integration Methods
% a
clc
clear

xu = 3.8;
xl = .8;
f = @(x) cot(x - (pi/4)) + exp(x);
method = 1;
n = 20;

Integ(f,xl,xu,method,n)

%%
% b
clc
clear

xu = 2*pi;
xl = 0;
f = @(x) sin(x) + cos(x);
method = 1;
n = 18;

Integ(f,xl,xu,method,n)

%% C RHS 
clc
clear

xu = 10;
xl = -10;
f = @(x) x.^2;
method = 2;
n = 20;

Integ(f,xl,xu,method,n)

%% C LHS
clc
clear

xu = 10;
xl = 1;
f = @(x) log(x);
method = 3;
n = 20;

Integ(f,xl,xu,method,n)

%% Integration Method Function

 
function Integ(f,xl,xu,method,n)
    % To use this function, input an anonymous function, lower bound, upper bound,integration method, and number of rectangles
    % Method = 1: Midpoint, 2: RHS, 3: LHS

    % The funciton will calculate the total area of all of the rectangles and
    % will output a plot to help visualize how the sum was computed
    
    int_ans = integral(f,xl,xu); % computes actual integral to later compare to calculated integral
    fplot(f,[xl-1,xu+1]) % this block plots the function and bounds of integration 
    hold on
    xline(xl);
    xline(xu);
    xlabel('x-axis')
    ylabel('y-axis')
 
    hold on
    total_area = 0;
    l = (xu-xl)/n; % width of rectangles
    span = linspace(xl,xu,n + 1); 
    if method == 1 % Midpoint Method calculation
        for ii = span(1:end-1) 
            hold on
            title('Midpoint Method')
            hold on 
            plot(ii+l/2,f(ii+l/2),'k.')
            hold on 
            xrec = [ii,ii,ii+l,ii+l];
            yrec = [0,f(ii+l/2),f(ii+l/2),0];
            a = area(xrec,yrec);
            if f(ii+l/2) < 0 % sums area of rectangles
                total_area = total_area - polyarea(xrec,yrec); % subtracts if negative
            else
                total_area = total_area + polyarea(xrec,yrec); % adds otherwise
            end
            set(a,'FaceColor','r');
        end

    elseif method == 2 % RHS Method calculation
        for ii = span(1:end-1)
            hold on
            title('Right Hand Sum Method')
            hold on
            plot(ii,f(ii),'k.')
            hold on 
            xrec = [ii,ii,ii+l,ii+l];
            yrec = [0,f(ii + l),f(ii + l),0];
            a = area(xrec,yrec);
            if f(ii + l) < 0 % sums area of rectangles. 
                total_area = total_area - polyarea(xrec,yrec); % subtracts if negative
            else
                total_area = total_area + polyarea(xrec,yrec); % adds otherwise
            end
            set(a,'FaceColor','r')
        end

    elseif method == 3 % LHS Method Calculation
        for ii = span(1:end-1)
            hold on
            title('Left Hand Sum Method')
            plot(ii,f(ii),'k.')
            hold on 
            xrec = [ii,ii,ii+l,ii+l];
            yrec = [0,f(ii),f(ii),0];
            a = area(xrec,yrec);
            if f(ii+l/2) < 0 % sums area of rectangles
                total_area = total_area - polyarea(xrec,yrec); % subtracts if negative
            else
                total_area = total_area + polyarea(xrec,yrec); % adds otherwise
            end
            set(a,'FaceColor','r');
        end
    else 
        disp('ERROR: Please input a valid integration method (1,2, or3).\n')
    end
    legend('function','lower bound','upper bound','points','Area Rectangles','location','eastoutside')
    fprintf('The numerical integral is %f units \n',int_ans)
    fprintf('The calculated area is %f units \n',total_area)
end