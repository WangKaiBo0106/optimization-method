close all;
clc;
N = 100;
D = 10;
T = 200;
c1 = 1.5;
c2 = 1.5;
Wmax = 0.8;
Wmin = 0.4;
w = 0.8;
Xmax = 20;
Xmin = -20;
Vmax = 20;
Vmin = -10;


%初始化种群个体
x = rand(N,D) * (Xmax - Xmin) + Xmin;	%随机生成N个*D维
v = rand(N,D) * (Vmax - Vmin) + Vmin;

%初始化个体最优位置和最优值
p = x;			%p--个体最优位置，N个
pbest = ones(N,1);	%pbest--个体最优值，N个
for i = 1:N
    pbest(i) = func1(x(i,:));
end
%初始化全局最优位置和最优值
g = ones(1,D);		%g--全局最优位置，1个
gbest = inf;		%gbest--全局最优值，1个
for i = 1:N
    if(pbest(i) < gbest)
        g = p(i,:);
	gbest = pbest(i);
    end
end
gb = ones(1,T);

%按照公式迭代
for i = 1:T
    for j = 1:N
	%更新个体最优位置和最优值
	if (func1(x(j,:)) < pbest(j))
	    p(j,:) = x(j,:);
	    pbest(j) = func1(x(j,:));
	end
	%更新全局最优位置和最优值
	if (pbest(j) < gbest)
	    g = p(j,:);
	    gbest = pbest(j,:);
	end

	%计算动态惯性权重
	%w = Wmax - (Wmax -Wmin)*i/T; %逐渐减小

	%更新位置和速度值
	v(j,:) = w*v(j,:) + c1*rand*(p(j,:)-x(j,:)) + c2*rand*(g-x(k,:));
	x(j,:) = x(j,:) + v(j,:);
	%边界条件处理
	for ii = 1:D
	    if (v(j,ii) > Vmax) || (v(j,ii) < Vmin)
		v(j,ii) = rand * (Vmax - Vmin) + Vmin;
	    end
	    if (x(j,ii) > Xmax) || (x(j,ii) < Xmin)
		x(j,ii) = rand * (Xmax - Xmin) + Xmin;
	    end
	end
    end
    %记录历代全局最优值
    gb(i) = gbest;
end

%g;
%gb(end);
figure
plot(gb)
xlabel('迭代次数');
ylabel('适应度值');
title('适应度进化曲线')

%适应度函数
function result = func1(x)
summ = sum(x.^2);
result = summ;
end
