close all;
clc;

C = 20;                             %总测试次数
T = 100;                            %每次迭代次数
D_list = [2 30 100];
N_list = [20 40 60];

% Sphere % Rosenbrock % Dropwave % Shubert % Rastrigin % Ackley
func = @Rosenbrock;
xMax = 10;      
xMin = -xMax;
target = 0;                         %全局最优值
allow = 100;                         %允许误差

select = 3;
D = D_list(select);                 %函数维度
N = N_list(select);                 %种群数目
beta = 0.2;                         %成年狮所占比例因子
step = 0.05 * (xMax - xMin);
localSearch = 10;                   %狮王局部搜索次数

gbtemp = zeros(C,T);                %每次测试、每次迭代的结果
xKingTrace = zeros(T,D);

%%%%%%%%%%%%%%%%%%%% 狮群状态初始化 %%%%%%%%%%%%%%%%%%%%
%初始化狮群中个体位置
xInicial = rand(N,D) * (xMax - xMin) + xMin;
%初始化 个体历史最优位置xBest 和 群体历史最优位置gBest
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
    for i = 1:T
        alphaF = step * exp(-30*((i/T)^10));    %母狮移动范围扰动因子
        alphaC = step * (T - i) / T;            %幼狮移动范围扰动因子
        

        %更新狮王位置
        for j = 1:localSearch
            alphaK = step * exp(-30*((j/localSearch)^10));
            
            xKing = gBest + (rand(1,D) * 2 - 1)*alphaK;
            xKingValue = func(xKing);
            if xKingValue < gBestValue
               gBest = xKing;
               gBestValue = xKingValue;
            end
        end
        %更新母狮位置
        for j = 1:n
            xHunter(j,:) = ((xBest(j,:)+xBest(ceil(rand*n),:))/2) * (1 + alphaF * rand);
            %边界条件处理
            for t = 1:D
                if (xHunter(j,t) > xMax) || (xHunter(j,t) < xMin)
                    xHunter(j,t) = rand * (xMax - xMin) + xMin;
                end
            end
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
            %边界条件处理
            for t = 1:D
                if (xCub(j,t) > xMax) || (xCub(j,t) < xMin)
                    xCub(j,t) = rand * (xMax - xMin) + xMin;
                end
            end
        end
        
         
        %%%%% 计算适应度,更新历史最优位置 %%%%%
        
        %更新个体历史最优位置xBest前半段-母狮
        xHunterValue = ones(n,1);	
        for t = 1:n
            xHunterValue(t) = func(xHunter(t,:));   
            if xHunterValue(t) < xBestValue(t)
                xBest(t,:) = xHunter(t,:);
                xBestValue(t) = xHunterValue(t);
            end
        end
        %更新个体历史最优位置xBest后半段-幼狮
        xCubValue = ones(N-n,1);
        for t = 1:N-n
            xCubValue(t) = func(xCub(t,:));
            if xCubValue(t) < xBestValue(n + t)
                xBest(n + t,:) = xCub(t,:);
                xBestValue(n + t) = xCubValue(t);
            end
        end       
        %更新母狮群体历史最优位置xHunterBest
        for t = 1:n
            if xHunterBestValue < xHunterValue(t)
                xHunterBest = xHunter(t,:);
                xHunterBestValue = xHunterValue(t);
            end
        end       
        %更新gBest
        for t = 1:N
            if xBestValue(t) < gBestValue
                gBest = xBest(t,:);
                gBestValue = xBestValue(t);
            end
        end
        
        gbtemp(k,i) = gBestValue;
        
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        % 抽样一次试验（迭代T次）记录狮群路径
        if k == C
            xKingTrace(i,:) = xKing;       
        end
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
    end
end

gb = mean(gbtemp);

%%%%%%%%%%%%%%%%%%%% 输出结果 %%%%%%%%%%%%%%%%%%%%
%狮群捕食路径(二维)
% figure
% plot(xKingTrace(:,1), xKingTrace(:,2),'ro-');
% figure
% for j = 1:n
%     plot(xBest(j,1), xBest(j,2),'ro');
%     hold on;
% end
% for j = n+1:N
%     plot(xBest(j,1), xBest(j,2),'ko');
%     hold on;
% end

%数据分析：成功率、收敛所需（平均）迭代次数、误差范围
sucRate = 0;
convTime = 0;
for i = 1:C
    if gbtemp(i,T) - target < allow     %达到迭代次数上限时满足要求则成功
        sucRate = sucRate + 1;
    end
end
fprintf('成功率：%.2f\n', sucRate / C);
fprintf('最优值：%e\n', min(gbtemp(:,T)));
fprintf('最差值：%e\n', max(gbtemp(:,T)));
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


%       问题
%   *1.局部搜索能力不够
%   *2.狮王的搜索结果未参与比较
%   *3.如何实现狮王在最优值附近的小范围搜索以提高收敛速度
%   4.母狮的随机性太大，和幼狮作用区分不大
%   *5.狮王局部搜索的步长
%   *6.位置更新后可能会超过限定范围，需要处理边界条件
%   7.绘制100次迭代过程中整体狮群位置的变化（二维）
%   *8.绘制100次迭代过程中狮王的移动路径
%   9.狮王的局部搜索应与总的迭代次数也有关系

%       解决方法
%   1.狮王在最优值附近小范围搜索，寻找更优点，和母狮幼狮的搜索结果比较
%     情况一：狮王发现了更优点 -- 及时更新gBest
%     情况二：母狮或者幼狮发现了更优点 -- 及时更新狮王位置
%   2.增加对狮王的配合或者狮王的小范围搜索频次，加快收敛






%%%%%%%%%%%%%%%%%%%% 适应度函数 %%%%%%%%%%%%%%%%%%%%
                  % 测试函数func %
function [y] = Sphere(xx)          %%%% Sphere函数，-10~10
    d = length(xx);
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + xi^2;
    end
    y = sum;
end

function [y] = Rosenbrock(xx)          %%%% Rosenbrock函数，-10~10
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

function [y] = Dropwave(xx)        	%%%% Dropwave函数，二元，-2~2
    x1 = xx(1);
    x2 = xx(2);
    
    frac1 = 1 + cos(12*sqrt(x1^2+x2^2));
    frac2 = 0.5*(x1^2+x2^2) + 2;

    y = -frac1/frac2;
end

function [y] = Shubert(xx)            %%%% Shubert函数，二元，-10~10
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


function [y] = Rastrigin(xx)            %%%% Rastrigin函数  -5.12~5.12
    d = length(xx);                 %min = 0
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + (xi^2 - 10*cos(2*pi*xi));
    end

    y = 10*d + sum;
end

function [y] = Ackley(xx, a, b, c)    %%%% Ackley函数  -32~32
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