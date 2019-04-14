function [sucRate,minv,maxv,meanv,stdv,convTime] = PSO(func,xMax,xMin,target,error,C,T,D,N,~,~,fitness)
    c1 = 2;
    c2 = 2;
    Wmax = 0.8;
    Wmin = 0.4;
%     w = 0.7298;

    Vmax = 20;
    Vmin = -10;
    gbtemp = zeros(C,T);                %ÿ�β��ԡ�ÿ�ε����Ľ��
    for k = 1:C
        %��ʼ����Ⱥ����
        x = rand(N,D) * (xMax - xMin) + xMin;	%�������N��*Dά
        v = rand(N,D) * (Vmax - Vmin) + Vmin;

        %��ʼ����������λ�ú�����ֵ
        p = x;			%p--��������λ�ã�N��
        pbest = ones(N,1);	%pbest--��������ֵ��N��
        for i = 1:N
            pbest(i) = func(x(i,:));
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

        %���չ�ʽ����
        for i = 1:T
            for j = 1:N
            %���¸�������λ�ú�����ֵ
            if (func(x(j,:)) < pbest(j))
                p(j,:) = x(j,:);
                pbest(j) = func(x(j,:));
            end
            %����ȫ������λ�ú�����ֵ
            if (pbest(j) < gbest)
                g = p(j,:);
                gbest = pbest(j,:);
            end

            %���㶯̬����Ȩ��
            w = Wmax - (Wmax -Wmin)*i/T; %�𽥼�С

            %����λ�ú��ٶ�ֵ
            v(j,:) = w*v(j,:) + c1*rand*(p(j,:)-x(j,:)) + c2*rand*(g-x(j,:));
            x(j,:) = x(j,:) + v(j,:);
            %�߽���������
            for ii = 1:D
                if (v(j,ii) > Vmax) || (v(j,ii) < Vmin)
                v(j,ii) = rand * (Vmax - Vmin) + Vmin;
                end
                if (x(j,ii) > xMax) || (x(j,ii) < xMin)
                x(j,ii) = rand * (xMax - xMin) + xMin;
                end
            end
            end
            %��¼����ȫ������ֵ
            gbtemp(k,i) = gbest;
        end
    end
    
    gb = mean(gbtemp);
    
    %���ݷ������ɹ��ʡ���Χ
    sucRate = 0;
    convTime = 0;
    for i = 1:C
        if gbtemp(i,T) - target < error     %�ﵽ������������ʱ����Ҫ����ɹ�
            sucRate = sucRate + 1;
        end
    end
    sucRate = sucRate / C;
    minv = min(gbtemp(:,T));
    maxv = max(gbtemp(:,T));
    meanv = mean(gbtemp(:,T));
    stdv = std(gbtemp(:,T));
    %�������裨ƽ������������
    for i = 1:T
        if gb(i) - target > error
            convTime = convTime + 1;
        end
    end  
    %��Ӧ�Ƚ�������
    if fitness == 1
        figure
        plot(gb)
        xlabel('��������');
        ylabel('��Ӧ��ֵf(x)');
        title(['��Ӧ�Ƚ�������,ά��D=',num2str(D)])
    end
    

end


