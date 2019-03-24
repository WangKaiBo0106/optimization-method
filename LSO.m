close all;
clc;

C = 20;                             %总测试次数
T = 100;                            %每次迭代次数

D = 10;                              %函数维度
xMax = 10;      
xMin = -10;
target = -186.7309;
allow = 1;

N = 20;                             %种群数目
beta = 0.2;                         %成年狮所占比例因子
step = 0.05 * (xMax - xMin);

gbtemp = zeros(C,T);                %所有迭代的结果

%%%%%%%%%%%%%%%%%%%% 狮群状态初始化 %%%%%%%%%%%%%%%%%%%%
%初始化狮群中个体位置
xInicial = rand(N,D) * (xMax - xMin) + xMin;
%初始化个体和群体最优位置
xBest = xInicial;                   %个体最优位置
xBestValue = ones(1,N);             %个体最优位置的适应度
for i = 1:N
    xBestValue(i) = func(xBest(i,:));
end
gBestValue = inf;       
for i = 1:N
    if xBestValue(i) < gBestValue
        gBestValue = xBestValue(i); %群体最优位置的适应度
        gBest = xBest(i,:);         %群体最优位置
    end
end

%%%%%%%%%%%%%%%%%%%% 狮群角色分配 %%%%%%%%%%%%%%%%%%%%
%狮王位置初始化
xKing = gBest;
xKingValue = gBestValue;
%母狮的位置初始化
n = floor(beta * N);
xHunter = xBest(1:n,:);
xHunterValue = xBestValue(1:n);
%母狮最优位置初始化
xHunterBestValue = inf;
for i = 1:n
    if xHunterBestValue > xHunterValue(i)
        xHunterBest = xHunter(i,:);
        xHunterBestValue = xHunterValue(i);
    end
end
%幼狮位置初始化
xCub = xBest(n+1:N,:);
xCubValue = xBestValue(n+1:N);

%%%%%%%%%%%%%%%%%%%% 迭代过程 %%%%%%%%%%%%%%%%%%%%
for k = 1:C
    iCounter = 0;
    for i = 1:T
        alphaF = step * exp(-30*((i/T)^10));    %母狮移动范围扰动因子
        alphaC = step * (T - i) / T;            %幼狮移动范围扰动因子

        %更新狮王位置
        xKing = gBest * (1 + rand * norm(xKing - gBest));
        %更新母狮位置
        for j = 1:n
            xHunter(j,:) = ((xBest(j,:)+xBest(ceil(rand*n),:))/2) * (1 + alphaF * rand);
        end
        %更新幼狮位置
        gWorst = xMin + xMax - gBest;	
        for j = 1:N-n
            q = rand;
            if q <= 1/3     %向狮王靠近
                xCub(j,:) = ((gBest + xBest(j + n,:))/2) * (1 + alphaC * rand);
            elseif q <2/3	%向母狮中的最优位置xHunterBest靠近
                xCub(j,:) = ((xHunterBest + xBest(j + n,:))/2) * (1 + alphaC * rand);
            else            %远离狮王
                xCub(j,:) = ((gWorst + xBest(j + n,:))/2) * (1 + alphaC * rand);
            end
        end

        %计算适应度，更新个体xBest 和 群体最优位置gBest
        % |-- 为保证收敛，只保留每个狮子的历史最优位置
        % |-- 如果位置更新后为更差的点，则放弃该位置
        
        %更新xBest前半段-母狮
        xHunterValue = ones(n,1);	
        for t = 1:n
            xHunterValue(t) = func(xHunter(t,:));   
            if xHunterValue(t) < xBestValue(t)
                xBest(t,:) = xHunter(t,:);
                xBestValue(t) = xHunterValue(t);
            end
        end
        %更新xBest后半段-幼狮
        xCubValue = ones(N-n,1);
        for t = 1:N-n
            xCubValue(t) = func(xCub(t,:));
            if xCubValue(t) < xBestValue(n + t)
                xBest(n + t,:) = xCub(t,:);
                xBestValue(n + t) = xCubValue(t);
            end
        end
        
        %更新gBest
        for t = 1:N
            if xBestValue(t) < gBestValue
                gBest = xBest(t,:);
                gBestValue = xBestValue(t);
            end
        end
        
        %更新xHunterBest
        for t = 1:n
            if xHunterBestValue < xHunterValue(t)
                xHunterBest = xHunter(t,:);
                xHunterBestValue = xHunterValue(t);
            end
        end
        
        %每迭代10次重新分配角色
        % |--10次之内，gBest每次都会更新，因此xKing只在初始阶段等于gBest
        iCounter = iCounter + 1;
        if iCounter == 10
            iCounter = 0;
            xKing = gBest;
        end
        gbtemp(k,i) = gBestValue;
    end
end

gb = mean(gbtemp);

%%%%%%%%%%%%%%%%%%%% 输出结果 %%%%%%%%%%%%%%%%%%%%
%数据分析
%成功率、收敛所需迭代次数、误差范围
sucRate = 0;
convTime = 0;
for i = 1:C
    if gbtemp(i,100) - target < allow
        sucRate = sucRate + 1;
    end
end
fprintf('成功率：%.2f\n', sucRate / C);
for i = 1:T
    if gb(i) - target > allow
        convTime = convTime + 1;
    end
end
fprintf('收敛所需迭代次数：%d\n', convTime);

%适应度进化曲线
figure
plot(gb)
xlabel('迭代次数');
ylabel('适应度值f(x)');
title(['适应度进化曲线,维度D=',num2str(D)])








%%%%%%%%%%%%%%%%%%%% 适应度函数 %%%%%%%%%%%%%%%%%%%%
         % 测试函数func %
function [y] = func1(xx)          %%%% Sphere函数，-10~10
    d = length(xx);
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + xi^2;
    end
    y = sum;
end

function [y] = func2(xx)          %%%% Rosenbrock函数，-10~10
    d = length(xx);
    sum = 0;
    for ii = 1:(d-1)
        xi = xx(ii);
        xnext = xx(ii+1);
        new = 100*(xnext-xi^2)^2 + (xi-1)^2;
        sum = sum + new;
    end
    y = sum;
end

function [y] = func3(xx)        	%%%% Dropwave函数，二元，-2~2
    x1 = xx(1);
    x2 = xx(2);
    
    frac1 = 1 + cos(12*sqrt(x1^2+x2^2));
    frac2 = 0.5*(x1^2+x2^2) + 2;

    y = -frac1/frac2;
end

function [y] = func(xx)            %%%% Shubert函数，二元，-10~10
    x1 = xx(1);                     %min = -186.7309
    x2 = xx(2);
    sum1 = 0;
    sum2 = 0;

    for ii = 1:5
        new1 = ii * cos((ii+1)*x1+ii);
        new2 = ii * cos((ii+1)*x2+ii);
        sum1 = sum1 + new1;
        sum2 = sum2 + new2;
    end

    y = sum1 * sum2;
end


function [y] = func5(xx)            %%%% Rastrigin函数  -5.12~5.12
    d = length(xx);                 %min = 0
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + (xi^2 - 10*cos(2*pi*xi));
    end

    y = 10*d + sum;
end

function [y] = func6(xx, a, b, c)    %%%% Ackley函数  -32~32
    d = length(xx);                  %min = 0

    if (nargin < 4)
        c = 2*pi;
    end
    if (nargin < 3)
        b = 0.2;
    end
    if (nargin < 2)
        a = 20;
    end

    sum1 = 0;
    sum2 = 0;
    for ii = 1:d
        xi = xx(ii);
        sum1 = sum1 + xi^2;
        sum2 = sum2 + cos(c*xi);
    end

    term1 = -a * exp(-b*sqrt(sum1/d));
    term2 = -exp(sum2/d);

    y = term1 + term2 + a + exp(1);
end