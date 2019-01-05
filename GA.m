%%%%%%%%%%%%��׼�Ŵ��㷨%%%%%%%%%%%%%%%%
close all;
clc;
NP = 50;
L = 20;
Pc = 0.8;
Pm = 0.1;
G = 100;
Xs = 10;    %����������
Xx = 0;     %����������
f = randi([0,1],NP,L);
nf = zeros(NP,L);
x = zeros(1,NP);
Fit = zeros(1,NP);

%�Ŵ��㷨ѭ��
for k = 1:G
	%�����ƽ���Ϊ��������ʮ����
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
	fBest = f(rr(1,1),:);		%���Ÿ��������
	xBest = x(rr(1,1));			%���Ÿ������״
	Fit = (Fit - minFit) / (maxFit - minFit);	%��һ��
	
	%ѡ�����̶ģ�                      *����*
	sum_Fit = sum(Fit);
	fitvalue = Fit./sum_Fit;	
	fitvalue = cumsum(fitvalue);        %��Ӧ������ۼӣ����ӳ�0~1���䣩
	ms = sort(rand(NP,1));              %����������ms����
	fiti = 1;
	newi = 1;
	while newi <= NP	
		if ms(newi) < fitvalue(fiti)    %����Ӧ�ȵĸ���ᱻ���ѡ��
			nf(newi,:) = f(fiti,:);     %ѡ����ɺ�Ⱥ���������ֲ���
			newi = newi + 1;
		else
			fiti = fiti + 1;
		end
	end

	%����
	for i = 1:2:NP
		p = rand;
		if p < Pc                   %���ݽ������ѡ��һ������Ⱦɫ��
			q = randi([0,1],1,L);   %���ѡ����н���Ļ��򣨶�㣩
			for j = 1:L             
				if q(j) == 1
					temp = nf(i+1,j);
					nf(i+1,j) = nf(i,j);
					nf(i,j) = temp;
				end
			end
		end
	end

	%����
	i = 1;
	while i <= round(NP * Pc)		%���ݱ������ȷ��Ⱥ���еı�����
		h = randi([1,NP],1,1);      %���ѡ����б���ĸ���
		for j = 1:round(L* Pc)		
			g = randi([1,L],1,1);	%���ѡ����б���Ļ���
			nf(h,g) =~ nf(h,g);
		end
		i = i + 1;
	end

	f = nf;
	f(1,:) = fBest;
	trace(k) = maxFit;	%��¼ÿһ����������Ӧ��
end
xBest;
figure
plot(trace)
xlabel('��������')
ylabel('Ŀ�꺯��ֵ')
title('��Ӧ�Ƚ�������')

%��Ӧ�Ⱥ���
function result = func1(x)
	fit = x + 10 * sin(5 * x) + 7 * cos(4 * x);
	result = fit;
end
