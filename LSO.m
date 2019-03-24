close all;
clc;

C = 20;                             %�ܲ��Դ���
T = 100;                            %ÿ�ε�������

D = 10;                              %����ά��
xMax = 10;      
xMin = -10;
target = -186.7309;
allow = 1;

N = 20;                             %��Ⱥ��Ŀ
beta = 0.2;                         %����ʨ��ռ��������
step = 0.05 * (xMax - xMin);

gbtemp = zeros(C,T);                %���е����Ľ��

%%%%%%%%%%%%%%%%%%%% ʨȺ״̬��ʼ�� %%%%%%%%%%%%%%%%%%%%
%��ʼ��ʨȺ�и���λ��
xInicial = rand(N,D) * (xMax - xMin) + xMin;
%��ʼ�������Ⱥ������λ��
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
    iCounter = 0;
    for i = 1:T
        alphaF = step * exp(-30*((i/T)^10));    %ĸʨ�ƶ���Χ�Ŷ�����
        alphaC = step * (T - i) / T;            %��ʨ�ƶ���Χ�Ŷ�����

        %����ʨ��λ��
        xKing = gBest * (1 + rand * norm(xKing - gBest));
        %����ĸʨλ��
        for j = 1:n
            xHunter(j,:) = ((xBest(j,:)+xBest(ceil(rand*n),:))/2) * (1 + alphaF * rand);
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
        end

        %������Ӧ�ȣ����¸���xBest �� Ⱥ������λ��gBest
        % |-- Ϊ��֤������ֻ����ÿ��ʨ�ӵ���ʷ����λ��
        % |-- ���λ�ø��º�Ϊ����ĵ㣬�������λ��
        
        %����xBestǰ���-ĸʨ
        xHunterValue = ones(n,1);	
        for t = 1:n
            xHunterValue(t) = func(xHunter(t,:));   
            if xHunterValue(t) < xBestValue(t)
                xBest(t,:) = xHunter(t,:);
                xBestValue(t) = xHunterValue(t);
            end
        end
        %����xBest����-��ʨ
        xCubValue = ones(N-n,1);
        for t = 1:N-n
            xCubValue(t) = func(xCub(t,:));
            if xCubValue(t) < xBestValue(n + t)
                xBest(n + t,:) = xCub(t,:);
                xBestValue(n + t) = xCubValue(t);
            end
        end
        
        %����gBest
        for t = 1:N
            if xBestValue(t) < gBestValue
                gBest = xBest(t,:);
                gBestValue = xBestValue(t);
            end
        end
        
        %����xHunterBest
        for t = 1:n
            if xHunterBestValue < xHunterValue(t)
                xHunterBest = xHunter(t,:);
                xHunterBestValue = xHunterValue(t);
            end
        end
        
        %ÿ����10�����·����ɫ
        % |--10��֮�ڣ�gBestÿ�ζ�����£����xKingֻ�ڳ�ʼ�׶ε���gBest
        iCounter = iCounter + 1;
        if iCounter == 10
            iCounter = 0;
            xKing = gBest;
        end
        gbtemp(k,i) = gBestValue;
    end
end

gb = mean(gbtemp);

%%%%%%%%%%%%%%%%%%%% ������ %%%%%%%%%%%%%%%%%%%%
%���ݷ���
%�ɹ��ʡ��������������������Χ
sucRate = 0;
convTime = 0;
for i = 1:C
    if gbtemp(i,100) - target < allow
        sucRate = sucRate + 1;
    end
end
fprintf('�ɹ��ʣ�%.2f\n', sucRate / C);
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








%%%%%%%%%%%%%%%%%%%% ��Ӧ�Ⱥ��� %%%%%%%%%%%%%%%%%%%%
         % ���Ժ���func %
function [y] = func1(xx)          %%%% Sphere������-10~10
    d = length(xx);
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + xi^2;
    end
    y = sum;
end

function [y] = func2(xx)          %%%% Rosenbrock������-10~10
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

function [y] = func3(xx)        	%%%% Dropwave��������Ԫ��-2~2
    x1 = xx(1);
    x2 = xx(2);
    
    frac1 = 1 + cos(12*sqrt(x1^2+x2^2));
    frac2 = 0.5*(x1^2+x2^2) + 2;

    y = -frac1/frac2;
end

function [y] = func(xx)            %%%% Shubert��������Ԫ��-10~10
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


function [y] = func5(xx)            %%%% Rastrigin����  -5.12~5.12
    d = length(xx);                 %min = 0
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + (xi^2 - 10*cos(2*pi*xi));
    end

    y = 10*d + sum;
end

function [y] = func6(xx, a, b, c)    %%%% Ackley����  -32~32
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