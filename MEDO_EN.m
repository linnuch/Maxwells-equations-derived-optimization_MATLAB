function MEDO_EN()
%Data:2020/09/03
%Read me:This is a sample matlab code for the Maxwell's Equations Derived
%Optimization (MEDO). 

% Contact Information:
% Email: lll_work@buaa.edu.cn

% DISCLAIMER: This sample code is provided for educational purposes only. 
% Please read the referred paper to better understand the usage and the
% inner-workings of the algorithm.
% Use at own your risk! There is no guarantee that the code is bug free.

% Please refer to the following article in your publications:
%SU Donglin, LI Lilin, YANG Shunchuan, LI Bing, CHEN Guangzhi, XU Hui, 
%A new optimization algorithm applied in electromagnetics ¡ª Maxwell's equations 
%derived optimization (MEDO), In Journal of SCIENCE CHINA Information Sciences, 
%Volume 63, Issue 10, 2020, Pages 200301-, ISSN 1674-733X, 
%https://doi.org/10.1007/s11432-020-2927-2.
%(https://engine.scichina.com/doi/10.1007/s11432-020-2927-2)

%% preset of the coefficients
%more information about the physicle meaning and the recommended range of the coefficients 
%please refer to the article above. 
B = 1.16e-9;
pS = 563;
R1 = 1000;
g = 0.002;
L2 = 3e9;

%% preset of the population
popsize = 40;    % population size.
npar = 40;		% Dimension of the problem.
param.maxit = 500;		% Maximum number of iterations.


%% the statement and initialization of the variables
x = zeros(npar,popsize); %statement of the population
fx = zeros(1,popsize);  %statement of the objective function value
fxmin_global = zeros(1,param.maxit);  %statement of the global optimal solution
fxmin_f = zeros(1,param.maxit);  %statement of the optimal solution in the current iteration
u = zeros(npar,popsize);  %statement of the value of the voltage sourve
i2 = zeros(npar,popsize); %statement of the current of the loop1
i3 = zeros(npar,popsize); %statement of the current of the loop2
deltax = 0.0001; %calculation accuracy of difference method

lb = -ones(npar,1) * 5.12; %Lower dimension boundary.
ub = ones(npar,1) * 5.12; %Upper dimension boundary.
for i = 1:npar
       x(i,:) = lb(i) + (ub(i) - lb(i)) * rand(1, popsize);
end   %Randomize population
x_o = 0; %population in the last iteration, which is initialized to be 0 at the beginning.


%% calculate the objective funcion of the initialized polulation
for i = 1:popsize
    fx(i) = benchmark(x(:,i),npar);
end   %benchmark can be replaced with your own cost function and make sure the mapping above can carry out properly.

[fmin,xminnum] = min(fx);  % finding the best solution in the initial population :
fxmin_f(1) = fmin; % record the best solution in the current iteration
xmin = x(:, xminnum); % finding the best variable in the current iteration
glomin = fmin; %the global optimal solution in the whole interations
fxmin_global(1) = glomin; %record the global optimal solution
xminout = xmin;  % record the best variable in the whole iterations


fxlb = benchmark(lb,npar); % calculate the objective function of the lower boundary.
fmax = max(max(abs(fx)),abs(fxlb));  % the maximum value of the objective funcion in the current iteraion
fx=fx/fmax;  % normalize the objective functions of the population.
fxlb=fxlb/fmax; % normalize the objective function of the lower boundary.
h1 = 1 - fxlb; % the variable to calculate S3
h2 = 1 - fx(1); % the variable to calculate S3
R2 = 0.2; % the resistant of loop1
S3 =  0.1*ones(npar,1); % the loop inductance of loop2 


%% Start iterations :
for iter = 2:param.maxit  
    dx = x - x_o;
    x_o = x;
    for i = 1 : popsize  
        S3 = S3 + (h1 + h2) .* (x(:,1) - lb) ./ 2;
        for j = 2:i 
                h1 = 1 - fx(j -1);
                h2 = 1 - fx(j);           
                S3 = S3 + (h1 + h2) .* (x(: , j) - x(: , j-1)) ./ 2;  
        end  %calculation the area of loop2: S3
        R3 = 1-fx(i) + sqrt((x(:,i) - lb).^2+(fx(i)-fxlb).^2) + 0.1;  %calculation R3
        % the value of the voltage source is set to be the gradient of the objective function
        u(:,i) =  -derivation(x(:,i) , deltax ,npar); 
        u_max=max(abs(u(:,i))); 
        u(:,i)=u(:,i)/u_max;  %normalize the voltage source
 
        if iter == 2 % in the first iteration, i3 is calculated in condition of static field
            Ztotal = R1 + x(:,i)-lb + R2 .* R3 ./ (R2 + R3);
            i3(:,i) = u(:,i) ./ Ztotal ./ (1 + R3 ./ R2);
            i2(:,i) = R3 ./ R2 .* i3(:,i);
            i3oo = i3(:,i); %record i3 
            i2oo = i2(:,i); %record i2
        else 
            % in the later iterations, i3 is calculated in condition of
            % time-varying field
            [i2n,i3n] = solve_i(L2,dx(:,i),1-fx(i),u(:,i),i2(:,i),i2oo,i3(:,i),i3oo,R1 + x(:,i),R2,R3,B,S3);
            i2oo = i2(:,i); %record i2
            i3oo = i3(:,i); %record i3
            i2(:,i) = i2n; %update i2
            i3(:,i) = i3n; %update i3
        end
        x(:,i) = x(:,i) + (i3(:,i) * B / pS - (x(:,i) - 0) * g) * iter;  %update the population
    end 

    x = max(x,lb); 
    x = min(x,ub); %check the population
    for i = 1:popsize
        fx(i) = benchmark(x(:,i),npar);
    end  %update the values of the objective functions
    [fmin,xminnum] = min(fx);  %finding the best solution in the current iteration
    xmin = x(:, xminnum);  %finding the best variable in the current iteration
    glomin = min(fmin, glomin); %update the global optimal solution in the whole interations
    if fmin == glomin
        xminout = xmin; % update the best variable in the whole iterations
    end 
    fxmin_global(iter) = glomin; %record the global optimal solution
    fxmin_f(iter) = fmin;  % record the best solution in the current iteration
    fmax=max(max(abs(fx)),abs(fxlb));
    fx=fx/fmax;  %normalization
end


%% result
best_ans = min(fxmin_global)  %print the optimal solution
plot(fxmin_global); 
xlabel('the number of iteration');
ylabel('optimal solution');
title('convergence curve of MEDO');
end
 

%% the calculation of the current i3 in condition of time-varying field
%more information please refer to the referrence.
function [i2n,i3n] = solve_i(L2,dx,h,dfx,i2o,i2oo,i3o,i3oo,R1,R2,R3,B,S3) 
deltaS = h * dx + 0.5 * dfx .* dx .^ 2; % deltaS is approximately calculated by the area of a trapezoid.
L3 = S3; 
Z2 = abs(R2 + L2 .* (log(i2o) - log(i2oo)));
% the value of abs(1/i2*(di2/dt)) is transformed to be d(abs(log(i2)))/dt, 
% and the gradient is approximately calculated through the difference.
Z3 = abs(R3 + L3 .* (log(i3o) - log(i3oo)));
Ztotal = R1 + Z2 .* Z3 ./ (Z2 + Z3);
i3n = (dfx ./ Ztotal .* Z2 + B * deltaS) ./ (Z2 + Z3);
i2n = dfx ./ Ztotal - i3n;
i3n = (i3n) ./ (abs(i3n) + abs(i2n)); %normalization
i2n = (i2n) ./ (abs(i3n) + abs(i2n));
% i3n = (i3n) ./ (max(abs(i3n),abs(i2n))); %normalization
% i2n = (i2n) ./ (max(abs(i3n),abs(i2n)));
% end
end


%% approximately calculate the gradient through the difference method
% notice that if the benchmark is replaced with your own cost function, the
% input parameters of this function should be mapped.
function df = derivation(x,deltax,npar)
    y_n= zeros(npar,1);
    for i = 1:npar   
        x_n=x;
        x_n(i) = x(i) + deltax;
        y_n(i) = benchmark(x_n,npar); 
    end
    y = benchmark(x,npar);
    df = (y_n - y) / (deltax);
end


%% test funcion, which can be replaced with your own cost function.
function y = benchmark(x,npar) 
    y = sum(x.^2 - 10 * cos(2* pi *x) + 10);
end