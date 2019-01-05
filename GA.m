%%%%%%%%%%%%标准遗传算法%%%%%%%%%%%%%%%%
close all;
clc;
NP = 50;
L = 20;
Pc = 0.8;
Pm = 0.1;
G = 100;
Xs = 10;    %定义域上限
Xx = 0;     %定义域下限
f = randi([0,1],NP,L);
nf = zeros(NP,L);
x = zeros(1,NP);
Fit = zeros(1,NP);

%遗传算法循环
for k = 1:G
	%二进制解码为定义域内十进制
	for i = 1:NP
		U = f(i,:);
		m = 0;
		for j = 1:L
			m = U(j) * 2^(j-1) + m;
		end
		x(i) = Xx + m * (Xs - Xx) / (2^L-1);
		Fit(i) = func1(x(i));
	end
	maxFit = max(Fit);
	minFit = min(Fit);
	rr = find(Fit == maxFit);
	fBest = f(rr(1,1),:);		%最优个体的序列
	xBest = x(rr(1,1));			%最优个体的性状
	Fit = (Fit - minFit) / (maxFit - minFit);	%归一化
	
	%选择（轮盘赌）                      *方法*
	sum_Fit = sum(Fit);
	fitvalue = Fit./sum_Fit;	
	fitvalue = cumsum(fitvalue);        %适应度逐个累加（连接成0~1区间）
	ms = sort(rand(NP,1));              %必须对随机点ms排序
	fiti = 1;
	newi = 1;
	while newi <= NP	
		if ms(newi) < fitvalue(fiti)    %高适应度的个体会被多次选择
			nf(newi,:) = f(fiti,:);     %选择完成后群体数量保持不变
			newi = newi + 1;
		else
			fiti = fiti + 1;
		end
	end

	%交叉
	for i = 1:2:NP
		p = rand;
		if p < Pc                   %根据交叉概率选择一对相邻染色体
			q = randi([0,1],1,L);   %随机选择进行交叉的基因（多点）
			for j = 1:L             
				if q(j) == 1
					temp = nf(i+1,j);
					nf(i+1,j) = nf(i,j);
					nf(i,j) = temp;
				end
			end
		end
	end

	%变异
	i = 1;
	while i <= round(NP * Pc)		%根据变异概率确定群体中的变异数
		h = randi([1,NP],1,1);      %随机选择进行变异的个体
		for j = 1:round(L* Pc)		
			g = randi([1,L],1,1);	%随机选择进行变异的基因
			nf(h,g) =~ nf(h,g);
		end
		i = i + 1;
	end

	f = nf;
	f(1,:) = fBest;
	trace(k) = maxFit;	%记录每一代的最优适应度
end
xBest;
figure
plot(trace)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')

%适应度函数
function result = func1(x)
	fit = x + 10 * sin(5 * x) + 7 * cos(4 * x);
	result = fit;
end

