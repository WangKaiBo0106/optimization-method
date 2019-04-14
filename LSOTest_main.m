% Sphere        单峰 [-10,10]       f(0,...,0)=0          
% Rosenbrock    谷形 [-10,10]       f(1,...,1)=0         
% Rastrigin     多峰 [-5.12,5.12]   f(0,...,0)=0          
% Ackley        多峰 [-32,32]       f(0,...,0)=0  
% Dropwave      多峰 [-2,2]         f(0,0)    =-1             
% Shubert       多峰 [-10,10]       f(*,*)    =-186.7309   

close all;
clc;

D_list = [2, 30, 100];
N_list = [20, 40, 60];
Max_list = [10, 10, 2, 10, 5.12, 32];
Target_list = [0, 0, -1, -186.7309, 0, 0];
Allow_list = [1, 100, 0.01, 36, 0.01, 0.01];
Func_list = {@Sphere,@Rosenbrock,@Dropwave,@Shubert,@Rastrigin,@Ackley};
Funcname = {'Sphere','Rosenbrock','Dropwave','Shubert','Rastrigin','Ackley'};

test = 0;
C = 20;                                 %总测试次数
T = 200;                                %每次迭代次数
trac = 0;                               %是否绘制狮群轨迹
fitness = 0;                            %是否绘制适应度进化曲线

if test == 1
    %%%%%%%%%%%%%%%%%%%% 单元测试 %%%%%%%%%%%%%%%%%%%%
    Fs = 6;
    Ds = 1;
    func = Func_list{Fs};
    xMax = Max_list(Fs);    
    xMin = -xMax;
%     xMax = 15;
%     xMin = -5;
    target = Target_list(Fs);       %全局最优值
    error = Allow_list(Fs);         %允许误差
    D = D_list(Ds);                 %函数维度
    N = N_list(Ds);                 %种群数目
    beta = 0.2;                     %成年狮所占比例因子

    [sucRate,minv,maxv,meanv,stdv,convTime] = AdaLSO_beta(func,xMax,xMin,...
        target,error,C,T,D,N,beta,trac,fitness);

    fprintf('测试函数：%s，维度：%d，种群规模：%d\n', Funcname{Fs}, D, N);
    fprintf('成功率：%.2f\n', sucRate);
    fprintf('最优值：%.2e\n', minv);
    fprintf('最差值：%.2e\n', maxv);
    fprintf('平均值：%.2e\n', meanv);
    fprintf('标准差：%.2e\n', stdv);
    fprintf('收敛所需迭代次数：%d\n', convTime);
else
    %%%%%%%%%%%%%%%%%%%% 完全测试 %%%%%%%%%%%%%%%%%%%%
    data = zeros(6,14);
    s = 1;
    for Fs = 1:6
        for Ds = 1:3
            func = Func_list{Fs};
            xMax = Max_list(Fs);   
            xMin = -xMax;
            target = Target_list(Fs);       %全局最优值
            error = Allow_list(Fs);         %允许误差
            D = D_list(Ds);                 %函数维度
            N = N_list(Ds);                 %种群数目
            beta = 0.2;                     %成年狮所占比例因子

            [sucRate,minv,maxv,meanv,stdv,convTime] = AdaLSO_beta(func,xMax,xMin,...
                target,error,C,T,D,N,beta,trac,fitness);

            data(1,s) = sucRate;
            data(2,s) = minv;
            data(3,s) = maxv;
            data(4,s) = meanv;
            data(5,s) = stdv;
            data(6,s) = convTime;
            s = s + 1;

    %         fprintf('测试函数：%s，维度：%d，种群规模：%d\n', Funcname{Fs}, D, N);
    %         fprintf('成功率：%.2f\n', sucRate);
    %         fprintf('最优值：%.2e\n', minv);
    %         fprintf('最差值：%.2e\n', maxv);
    %         fprintf('平均值：%.2e\n', meanv);
    %         fprintf('标准差：%.2e\n', stdv);
    %         fprintf('收敛所需迭代次数：%d\n', convTime);
    %         fprintf('--------------------------------------------\n');

            if Fs == 3 || Fs == 4
                break;
            end
        end
    end
    xlswrite('data.xlsx',data,1,'B37:O42');
end






%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%       问题
% *1.局部搜索能力不够
% *2.狮王的搜索结果未参与比较
% *3.如何实现狮王在最优值附近的小范围搜索以提高收敛速度
% *4.母狮的随机性太大，和幼狮作用区分不大
% *5.狮王局部搜索的步长
% *6.位置更新后可能会超过限定范围，需要处理边界条件
% *7.绘制迭代过程中整体狮群位置的变化（二维）
% *8.绘制100次迭代过程中狮王的移动路径
% *9.狮王的局部搜索应与总的迭代次数也有关系
% 10.制定基于种群分布熵的策略，指导行为：勘探or开发
% 11.母狮的全局搜索能力有限，依赖于初始值，缺乏勘探能力
% 12.当前的位置更新策略太过贪婪、短视，不考虑长期收益
% 13.当最优位置不在寻优空间中心时，算法性能急剧下降
%
%       解决方法
% 1.狮王在最优值附近小范围搜索，寻找更优点，和母狮幼狮的搜索结果比较
%     情况一：狮王发现了更优点 -- 及时更新gBest
%     情况二：母狮或者幼狮发现了更优点 -- 及时更新狮王位置
% 2.增加对狮王的配合或者狮王的小范围搜索频次，加快收敛
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%%%%%%%%%%%%%%%%%%%% 适应度函数 %%%%%%%%%%%%%%%%%%%%
function [y] = Sphere(xx)  
    d = length(xx);
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + xi^2;
    end
    y = sum;
end


function [y] = Rosenbrock(xx)   
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


function [y] = Dropwave(xx)      
    x1 = xx(1);
    x2 = xx(2);
    
    frac1 = 1 + cos(12*sqrt(x1^2+x2^2));
    frac2 = 0.5*(x1^2+x2^2) + 2;

    y = -frac1/frac2;
end


function [y] = Shubert(xx)       
    x1 = xx(1);                
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



function [y] = Rastrigin(xx)          
    d = length(xx);                 
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + (xi^2 - 10*cos(2*pi*xi));
    end

    y = 10*d + sum;
end


function [y] = Ackley(xx, a, b, c)    
    d = length(xx);                  

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


function [y] = griewank(xx)
d = length(xx);
sum = 0;
prod = 1;

for ii = 1:d
	xi = xx(ii);
	sum = sum + xi^2/4000;
	prod = prod * cos(xi/sqrt(ii));
end

y = sum - prod + 1;

end