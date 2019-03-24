close all;
clc;

C = 20;                             %�ܲ��Դ���
T = 100;                            %ÿ�ε�������
D_list = [2 30 100];
N_list = [20 40 60];

% Sphere % Rosenbrock % Dropwave % Shubert % Rastrigin % Ackley
func = @Rosenbrock;
xMax = 10;      
xMin = -xMax;
target = 0;                         %ȫ������ֵ
allow = 100;                         %�������

select = 3;
D = D_list(select);                 %����ά��
N = N_list(select);                 %��Ⱥ��Ŀ
beta = 0.2;                         %����ʨ��ռ��������
step = 0.05 * (xMax - xMin);
localSearch = 10;                   %ʨ���ֲ���������

gbtemp = zeros(C,T);                %ÿ�β��ԡ�ÿ�ε����Ľ��
xKingTrace = zeros(T,D);

%%%%%%%%%%%%%%%%%%%% ʨȺ״̬��ʼ�� %%%%%%%%%%%%%%%%%%%%
%��ʼ��ʨȺ�и���λ��
xInicial = rand(N,D) * (xMax - xMin) + xMin;
%��ʼ�� ������ʷ����λ��xBest �� Ⱥ����ʷ����λ��gBest
xBest = xInicial;                   %��������λ��
xBestValue = ones(1,N);             %��������λ�õ���Ӧ��
for i = 1:N
    xBestValue(i) = func(xBest(i,:));
end
gBestValue = inf;       
for i = 1:N
    if xBestValue(i) < gBestValue
        gBestValue = xBestValue(i); %Ⱥ������λ�õ���Ӧ��
        gBest = xBest(i,:);         %Ⱥ������λ��
    end
end

%%%%%%%%%%%%%%%%%%%% ʨȺ��ɫ���� %%%%%%%%%%%%%%%%%%%%
%ʨ��λ�ó�ʼ��
xKing = gBest;
xKingValue = gBestValue;
%ĸʨ��λ�ó�ʼ��
n = floor(beta * N);
xHunter = xBest(1:n,:);
xHunterValue = xBestValue(1:n);
%ĸʨ����λ�ó�ʼ��
xHunterBestValue = inf;
for i = 1:n
    if xHunterBestValue > xHunterValue(i)
        xHunterBest = xHunter(i,:);
        xHunterBestValue = xHunterValue(i);
    end
end
%��ʨλ�ó�ʼ��
xCub = xBest(n+1:N,:);
xCubValue = xBestValue(n+1:N);

%%%%%%%%%%%%%%%%%%%% �������� %%%%%%%%%%%%%%%%%%%%
for k = 1:C
    for i = 1:T
        alphaF = step * exp(-30*((i/T)^10));    %ĸʨ�ƶ���Χ�Ŷ�����
        alphaC = step * (T - i) / T;            %��ʨ�ƶ���Χ�Ŷ�����
        

        %����ʨ��λ��
        for j = 1:localSearch
            alphaK = step * exp(-30*((j/localSearch)^10));
            
            xKing = gBest + (rand(1,D) * 2 - 1)*alphaK;
            xKingValue = func(xKing);
            if xKingValue < gBestValue
               gBest = xKing;
               gBestValue = xKingValue;
            end
        end
        %����ĸʨλ��
        for j = 1:n
            xHunter(j,:) = ((xBest(j,:)+xBest(ceil(rand*n),:))/2) * (1 + alphaF * rand);
            %�߽���������
            for t = 1:D
                if (xHunter(j,t) > xMax) || (xHunter(j,t) < xMin)
                    xHunter(j,t) = rand * (xMax - xMin) + xMin;
                end
            end
        end
        %������ʨλ��
        gWorst = xMin + xMax - gBest;	
        for j = 1:N-n
            q = rand;
            if q <= 1/3     %��ʨ������
                xCub(j,:) = ((gBest + xBest(j + n,:))/2) * (1 + alphaC * rand);
            elseif q <2/3	%��ĸʨ�е�����λ��xHunterBest����
                xCub(j,:) = ((xHunterBest + xBest(j + n,:))/2) * (1 + alphaC * rand);
            else            %Զ��ʨ��
                xCub(j,:) = ((gWorst + xBest(j + n,:))/2) * (1 + alphaC * rand);
            end
            %�߽���������
            for t = 1:D
                if (xCub(j,t) > xMax) || (xCub(j,t) < xMin)
                    xCub(j,t) = rand * (xMax - xMin) + xMin;
                end
            end
        end
        
         
        %%%%% ������Ӧ��,������ʷ����λ�� %%%%%
        
        %���¸�����ʷ����λ��xBestǰ���-ĸʨ
        xHunterValue = ones(n,1);	
        for t = 1:n
            xHunterValue(t) = func(xHunter(t,:));   
            if xHunterValue(t) < xBestValue(t)
                xBest(t,:) = xHunter(t,:);
                xBestValue(t) = xHunterValue(t);
            end
        end
        %���¸�����ʷ����λ��xBest����-��ʨ
        xCubValue = ones(N-n,1);
        for t = 1:N-n
            xCubValue(t) = func(xCub(t,:));
            if xCubValue(t) < xBestValue(n + t)
                xBest(n + t,:) = xCub(t,:);
                xBestValue(n + t) = xCubValue(t);
            end
        end       
        %����ĸʨȺ����ʷ����λ��xHunterBest
        for t = 1:n
            if xHunterBestValue < xHunterValue(t)
                xHunterBest = xHunter(t,:);
                xHunterBestValue = xHunterValue(t);
            end
        end       
        %����gBest
        for t = 1:N
            if xBestValue(t) < gBestValue
                gBest = xBest(t,:);
                gBestValue = xBestValue(t);
            end
        end
        
        gbtemp(k,i) = gBestValue;
        
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        % ����һ�����飨����T�Σ���¼ʨȺ·��
        if k == C
            xKingTrace(i,:) = xKing;       
        end
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
    end
end

gb = mean(gbtemp);

%%%%%%%%%%%%%%%%%%%% ������ %%%%%%%%%%%%%%%%%%%%
%ʨȺ��ʳ·��(��ά)
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

%���ݷ������ɹ��ʡ��������裨ƽ����������������Χ
sucRate = 0;
convTime = 0;
for i = 1:C
    if gbtemp(i,T) - target < allow     %�ﵽ������������ʱ����Ҫ����ɹ�
        sucRate = sucRate + 1;
    end
end
fprintf('�ɹ��ʣ�%.2f\n', sucRate / C);
fprintf('����ֵ��%e\n', min(gbtemp(:,T)));
fprintf('���ֵ��%e\n', max(gbtemp(:,T)));
for i = 1:T
    if gb(i) - target > allow
        convTime = convTime + 1;
    end
end
fprintf('�����������������%d\n', convTime);


%��Ӧ�Ƚ�������
figure
plot(gb)
xlabel('��������');
ylabel('��Ӧ��ֵf(x)');
title(['��Ӧ�Ƚ�������,ά��D=',num2str(D)])


%       ����
%   *1.�ֲ�������������
%   *2.ʨ�����������δ����Ƚ�
%   *3.���ʵ��ʨ��������ֵ������С��Χ��������������ٶ�
%   4.ĸʨ�������̫�󣬺���ʨ�������ֲ���
%   *5.ʨ���ֲ������Ĳ���
%   *6.λ�ø��º���ܻᳬ���޶���Χ����Ҫ����߽�����
%   7.����100�ε�������������ʨȺλ�õı仯����ά��
%   *8.����100�ε���������ʨ�����ƶ�·��
%   9.ʨ���ľֲ�����Ӧ���ܵĵ�������Ҳ�й�ϵ

%       �������
%   1.ʨ��������ֵ����С��Χ������Ѱ�Ҹ��ŵ㣬��ĸʨ��ʨ����������Ƚ�
%     ���һ��ʨ�������˸��ŵ� -- ��ʱ����gBest
%     �������ĸʨ������ʨ�����˸��ŵ� -- ��ʱ����ʨ��λ��
%   2.���Ӷ�ʨ������ϻ���ʨ����С��Χ����Ƶ�Σ��ӿ�����






%%%%%%%%%%%%%%%%%%%% ��Ӧ�Ⱥ��� %%%%%%%%%%%%%%%%%%%%
                  % ���Ժ���func %
function [y] = Sphere(xx)          %%%% Sphere������-10~10
    d = length(xx);
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + xi^2;
    end
    y = sum;
end

function [y] = Rosenbrock(xx)          %%%% Rosenbrock������-10~10
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

function [y] = Dropwave(xx)        	%%%% Dropwave��������Ԫ��-2~2
    x1 = xx(1);
    x2 = xx(2);
    
    frac1 = 1 + cos(12*sqrt(x1^2+x2^2));
    frac2 = 0.5*(x1^2+x2^2) + 2;

    y = -frac1/frac2;
end

function [y] = Shubert(xx)            %%%% Shubert��������Ԫ��-10~10
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


function [y] = Rastrigin(xx)            %%%% Rastrigin����  -5.12~5.12
    d = length(xx);                 %min = 0
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + (xi^2 - 10*cos(2*pi*xi));
    end

    y = 10*d + sum;
end

function [y] = Ackley(xx, a, b, c)    %%%% Ackley����  -32~32
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