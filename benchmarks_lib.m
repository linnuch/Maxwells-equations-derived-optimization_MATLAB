function y = benchmarks_lib(num,x,npar)
%Data:2020/09/03
%Read me:This is some benchmarks that can be used to test the optimization
%algorithms.

%% unimodal benchmark functions
if num == 1  %[-100,100]
    y = sum((x-2).^2);
end

if num == 2     %[-10,10]
    temp = 1;
    for ii = 1:npar
        temp = temp * abs(x(ii));
    end
    y = sum(abs(x)) + temp;
end

if num == 3     %[-100,100]
    y = 0;
    for ii = 1: npar
        for jj = 1:ii
            y1 = sum(x(jj));
        end
        y1 = y1.^2;
        y = y + y1;
    end
end

if num == 4   %[-100,100]
    y = sum((x+0.5).^2);
end

if num == 5  %[-1.28,1.28]
    t1 = 0;
    for ii = 1:npar
        t1 = t1 + ii * x(ii)^4;
    end
    y = t1 + rand;
end


%% multimodal benchmark functions
if num ==  6  %0
    y = sum(x.^2 - 10 * cos(2* pi *x) + 10);
end

if num == 7   %0
    y = -20 * exp(-0.2 * sqrt(sum(x.^2)/npar)) - exp(sum(cos(2*pi*x))/npar) + 20 + exp(1);
end

if num ==8  %0
    y1 = 1;
    for ii = 1:npar
        y1 = y1 * cos(x(ii)/sqrt(ii));
    end
    y = sum(x.^2) / 4000 - y1 + 1;
end

%% fixed-dimension multimodal benchmark functions
if num==9  %3
    y1 = 1 + (x(1) + x(2) + 1)^2 * (19 - 14*x(1) + 3 * x(1)^2 - 14*x(2) + 6 * x(1) * x(2) + 3*x(2)^2);
    y2 = 30 + (2*x(1) - 3*x(2))^2 * (18 - 32*x(1) + 12*x(1)^2 + 48*x(2) - 36*x(1)*x(2) + 27 * x(2)^2);
    y = y1*y2;
end

if (num==10)%Schaffer function N. 2   0
    y = 0.5 + ((sin(x(1)^2-x(2)^2))^2 - 0.5) / (1+0.001*(x(1)^2+x(2)^2))^2;
end

if (num == 11) %Schaffer  -1
    x1=x(1);  
    y1=x(2);  
    temp=x1^2+y1^2;  
    f=0.5-(sin(sqrt(temp))^2-0.5)/(1+0.001*temp)^2;  
    y=-f;
end


end


