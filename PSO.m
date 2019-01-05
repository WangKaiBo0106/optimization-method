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


%��ʼ����Ⱥ����
x = rand(N,D) * (Xmax - Xmin) + Xmin;	%�������N��*Dά
v = rand(N,D) * (Vmax - Vmin) + Vmin;

%��ʼ����������λ�ú�����ֵ
p = x;			%p--��������λ�ã�N��
pbest = ones(N,1);	%pbest--��������ֵ��N��
for i = 1:N
    pbest(i) = func1(x(i,:));
end
%��ʼ��ȫ������λ�ú�����ֵ
g = ones(1,D);		%g--ȫ������λ�ã�1��
gbest = inf;		%gbest--ȫ������ֵ��1��
for i = 1:N
    if(pbest(i) < gbest)
        g = p(i,:);
	gbest = pbest(i);
    end
end
gb = ones(1,T);

%���չ�ʽ����
for i = 1:T
    for j = 1:N
	%���¸�������λ�ú�����ֵ
	if (func1(x(j,:)) < pbest(j))
	    p(j,:) = x(j,:);
	    pbest(j) = func1(x(j,:));
	end
	%����ȫ������λ�ú�����ֵ
	if (pbest(j) < gbest)
	    g = p(j,:);
	    gbest = pbest(j,:);
	end

	%���㶯̬����Ȩ��
	%w = Wmax - (Wmax -Wmin)*i/T; %�𽥼�С

	%����λ�ú��ٶ�ֵ
	v(j,:) = w*v(j,:) + c1*rand*(p(j,:)-x(j,:)) + c2*rand*(g-x(k,:));
	x(j,:) = x(j,:) + v(j,:);
	%�߽���������
	for ii = 1:D
	    if (v(j,ii) > Vmax) || (v(j,ii) < Vmin)
		v(j,ii) = rand * (Vmax - Vmin) + Vmin;
	    end
	    if (x(j,ii) > Xmax) || (x(j,ii) < Xmin)
		x(j,ii) = rand * (Xmax - Xmin) + Xmin;
	    end
	end
    end
    %��¼����ȫ������ֵ
    gb(i) = gbest;
end

%g;
%gb(end);
figure
plot(gb)
xlabel('��������');
ylabel('��Ӧ��ֵ');
title('��Ӧ�Ƚ�������')

%��Ӧ�Ⱥ���
function result = func1(x)
summ = sum(x.^2);
result = summ;
end
