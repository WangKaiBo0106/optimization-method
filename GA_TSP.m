%%%%%%%%%%遗传算法解决旅行商问题%%%%%%%%%%
close all;
clc;
C = [1304 2312;3639 1315;4177 2244;3712 1399;3488 1535;3326 1556;...
	 3238 1229;4196 1044;4312  790;4386  570;3007 1970;2562 1756;...
	 2788 1491;2381 1676;1332  695;3715 1678;3918 2179;4061 2370;...
	 3780 2212;3676 2578;4029 2838;4263 2931;3429 1908;3507 2376;...
	 3394 2643;3439 3201;2935 3240;3140 3550;2545 2357;2778 2826;...
	 2370 2975];
N = size(C,1);  %城市数目
D = zeros(N);	%任意两个城市之间的间隔矩阵

%%%%%%%%%%求任意两个城市之间的间隔矩阵%%%%%%%%%%
for i = 1:N
	for j = 1:N
		D(i,j) = ((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)^0.5;
	end
end

NP = 200;
G = 1000;
f = zeros(NP,N);
F = [];
for i = 1:NP
	f(i,:) = randperm(N);	%初试种群-生成1-N的排列，代表路径
end
R = f(1,:);                 %最优个体-最短路径
len = zeros(NP,1);          %记录个体的路径长度
fitness = zeros(NP,1);      %归一化适应度
Rlength = [1,G];            %记录群体最短路径长度变化
gen = 0;

while gen < G
	%%%%%%%%%%计算路径长度%%%%%%%%%%
	for i = 1:NP
		len(i) = D(f(i,N),f(i,1)); %回到起点的路径长度
		for j = 1:N-1
			len(i) = len(i) + D(f(i,j),f(i,j+1));
		end
	end
	maxlen = max(len);
	minlen = min(len);
	%%%%%%%%%%更新最优个体（路径）%%%%%%%%%%
	rr = find(len == minlen);
	R = f(rr(1,1),:);
	%%%%%%%%%%归一化适应度（路径长度越短 适应度越高）%%%%%%%%%%
	for i = 1:length(len)
		fitness(i) = 1 - ((len(i) - minlen) / (maxlen - minlen + 0.001)); 
	end


	%%%%%%%%%%选择%%%%%%%%%%
	%非轮盘赌
	%选择算子只要体现“适应度大的存活率更高即可”
	nn = 0;
	for i = 1:NP
		if fitness(i) >=rand 		%选择之后种群数量可能减少
			nn = nn + 1;
			F(nn,:) = f(i,:);
		end
	end

	%当选择使种群数量减少时
	%通过交叉变异产生新个体补充
	%保持种群个数不变
	[aa,bb] = size(F);	
	while aa < NP				
		nnper = randperm(nn);		
		A = F(nnper(1),:);			
		B = F(nnper(2),:);
        
		%%%%%%%%%%交叉%%%%%%%%%%
		%需要保持路径组合有效性
		%即：每个城市只能经过一次
		%先对应位置交叉，然后调整
		W = ceil(N / 10);
		p = unidrnd(N - W + 1);
		for i = 1:W
			x = find(A == B(1,p + i - 1));
			y = find(B == A(1,p + i - 1));
			temp = A(1,p + i -1);
			A(1,p + i -1) = B(1,p + i -1);
			B(1,p + i -1) = temp;
			temp = A(1,x);
			A(1,x) = B(1,y);
			B(1,y) = temp;
        end
        
		%%%%%%%%%%变异%%%%%%%%%%
		%也需要保持路径组合的有效性
		%通过交换个体内的一对城市坐标来代替任意方向的变异
		p1 = floor(1 + N*rand);
		p2 = floor(1 + N*rand);
		while p1 == p2
			p1 = floor(1 + N*rand);
			p2 = floor(1 + N*rand);
		end
		tmp = A(p1);
		A(p1) = A(p2);
		A(p2) = tmp;
		tmp = B(p1);
		B(p1) = B(p2);
		B(p2) = tmp;
		F = [F;A;B];    %将交叉变异个体补充到种群中
		[aa,bb] = size(F);
    end
    
    
	if aa > NP
		F = F(1:NP,:);	%剔除多余的个体，保持种群规模不变
	end
	f = F;
	f(1,:) = R;
	clear F;
	gen = gen + 1;
	Rlength(gen) = minlen;
end

figure
for i = 1:N-1
	plot([C(R(i),1),C(R(i+1),1)],[C(R(i),2),C(R(i+1),2)],'bo-');
	hold on;
end
plot([C(R(N),1),C(R(1),1)],[C(R(N),2),C(R(1),2)],'ro-');
title(['优化最短距离：',num2str(minlen)]);
figure
plot(Rlength)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')

