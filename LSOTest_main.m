% Sphere        ���� [-10,10]       f(0,...,0)=0          
% Rosenbrock    ���� [-10,10]       f(1,...,1)=0         
% Rastrigin     ��� [-5.12,5.12]   f(0,...,0)=0          
% Ackley        ��� [-32,32]       f(0,...,0)=0  
% Dropwave      ��� [-2,2]         f(0,0)    =-1             
% Shubert       ��� [-10,10]       f(*,*)    =-186.7309   

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
C = 20;                                 %�ܲ��Դ���
T = 200;                                %ÿ�ε�������
trac = 0;                               %�Ƿ����ʨȺ�켣
fitness = 0;                            %�Ƿ������Ӧ�Ƚ�������

if test == 1
    %%%%%%%%%%%%%%%%%%%% ��Ԫ���� %%%%%%%%%%%%%%%%%%%%
    Fs = 6;
    Ds = 1;
    func = Func_list{Fs};
    xMax = Max_list(Fs);    
    xMin = -xMax;
%     xMax = 15;
%     xMin = -5;
    target = Target_list(Fs);       %ȫ������ֵ
    error = Allow_list(Fs);         %�������
    D = D_list(Ds);                 %����ά��
    N = N_list(Ds);                 %��Ⱥ��Ŀ
    beta = 0.2;                     %����ʨ��ռ��������

    [sucRate,minv,maxv,meanv,stdv,convTime] = AdaLSO_beta(func,xMax,xMin,...
        target,error,C,T,D,N,beta,trac,fitness);

    fprintf('���Ժ�����%s��ά�ȣ�%d����Ⱥ��ģ��%d\n', Funcname{Fs}, D, N);
    fprintf('�ɹ��ʣ�%.2f\n', sucRate);
    fprintf('����ֵ��%.2e\n', minv);
    fprintf('���ֵ��%.2e\n', maxv);
    fprintf('ƽ��ֵ��%.2e\n', meanv);
    fprintf('��׼�%.2e\n', stdv);
    fprintf('�����������������%d\n', convTime);
else
    %%%%%%%%%%%%%%%%%%%% ��ȫ���� %%%%%%%%%%%%%%%%%%%%
    data = zeros(6,14);
    s = 1;
    for Fs = 1:6
        for Ds = 1:3
            func = Func_list{Fs};
            xMax = Max_list(Fs);   
            xMin = -xMax;
            target = Target_list(Fs);       %ȫ������ֵ
            error = Allow_list(Fs);         %�������
            D = D_list(Ds);                 %����ά��
            N = N_list(Ds);                 %��Ⱥ��Ŀ
            beta = 0.2;                     %����ʨ��ռ��������

            [sucRate,minv,maxv,meanv,stdv,convTime] = AdaLSO_beta(func,xMax,xMin,...
                target,error,C,T,D,N,beta,trac,fitness);

            data(1,s) = sucRate;
            data(2,s) = minv;
            data(3,s) = maxv;
            data(4,s) = meanv;
            data(5,s) = stdv;
            data(6,s) = convTime;
            s = s + 1;

    %         fprintf('���Ժ�����%s��ά�ȣ�%d����Ⱥ��ģ��%d\n', Funcname{Fs}, D, N);
    %         fprintf('�ɹ��ʣ�%.2f\n', sucRate);
    %         fprintf('����ֵ��%.2e\n', minv);
    %         fprintf('���ֵ��%.2e\n', maxv);
    %         fprintf('ƽ��ֵ��%.2e\n', meanv);
    %         fprintf('��׼�%.2e\n', stdv);
    %         fprintf('�����������������%d\n', convTime);
    %         fprintf('--------------------------------------------\n');

            if Fs == 3 || Fs == 4
                break;
            end
        end
    end
    xlswrite('data.xlsx',data,1,'B37:O42');
end






%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%       ����
% *1.�ֲ�������������
% *2.ʨ�����������δ����Ƚ�
% *3.���ʵ��ʨ��������ֵ������С��Χ��������������ٶ�
% *4.ĸʨ�������̫�󣬺���ʨ�������ֲ���
% *5.ʨ���ֲ������Ĳ���
% *6.λ�ø��º���ܻᳬ���޶���Χ����Ҫ����߽�����
% *7.���Ƶ�������������ʨȺλ�õı仯����ά��
% *8.����100�ε���������ʨ�����ƶ�·��
% *9.ʨ���ľֲ�����Ӧ���ܵĵ�������Ҳ�й�ϵ
% 10.�ƶ�������Ⱥ�ֲ��صĲ��ԣ�ָ����Ϊ����̽or����
% 11.ĸʨ��ȫ�������������ޣ������ڳ�ʼֵ��ȱ����̽����
% 12.��ǰ��λ�ø��²���̫��̰�������ӣ������ǳ�������
% 13.������λ�ò���Ѱ�ſռ�����ʱ���㷨���ܼ����½�
%
%       �������
% 1.ʨ��������ֵ����С��Χ������Ѱ�Ҹ��ŵ㣬��ĸʨ��ʨ����������Ƚ�
%     ���һ��ʨ�������˸��ŵ� -- ��ʱ����gBest
%     �������ĸʨ������ʨ�����˸��ŵ� -- ��ʱ����ʨ��λ��
% 2.���Ӷ�ʨ������ϻ���ʨ����С��Χ����Ƶ�Σ��ӿ�����
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


%%%%%%%%%%%%%%%%%%%% ��Ӧ�Ⱥ��� %%%%%%%%%%%%%%%%%%%%
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