function [sucRate,minv,maxv,meanv,stdv,convTime] = LSO(func,xMax,xMin,target,error,C,T,D,N,beta,trac,fitness)
    %�������
    %func:���Ժ���
    %xMax:����������
    %target:Ѱ��Ŀ��
    %error:�������
    %C:�������
    %T:����������
    %D:����ά��
    %N:��Ⱥ����
    %beta:����ʨ��ռ��������
    %trac:�Ƿ���ƶ�ά��Ⱥ·��
    %fitness:�Ƿ������Ӧ�Ƚ�������
    
    %����ֵ
    %sucRate:�ɹ���
    %min:Ѱ�Ž����Сֵ
    %max:Ѱ�Ž�����ֵ
    %mean:Ѱ�Ž��ƽ��ֵ
    %std:Ѱ�Ž����׼��
    %convTime:���������������
    
    step = 0.05 * (xMax - xMin);

    gbtemp = zeros(C,T);                %ÿ�β��ԡ�ÿ�ε����Ľ��
    xKingTrace = zeros(T,D);
    densTrac =zeros(1,T);

    if trac == 1 && D == 2
        figure
        s = 1;  %��ͼ
    elseif trac == 2 && D == 2
        figure
        aviObj = VideoWriter('AdaLSO.avi');%����Ϊavi
        aviObj.FrameRate = 20;
        open(aviObj);
    end
    %%%%%%%%%%%%%%%%%%%% ��ʼ���� %%%%%%%%%%%%%%%%%%%%
    for k = 1:C

        %%%%%%%%%%%%%%%%%%%% ʨȺ״̬��ʼ�� %%%%%%%%%%%%%%%%%%%%
        %��ʼ�� ������ʷ����λ��xBest �� Ⱥ����ʷ����λ��gBest
        xInicial = rand(N,D) * (xMax - xMin) + xMin;
        xBest = xInicial;                   
        xBestValue = ones(1,N);             %��������λ�õ���Ӧ��
        for i = 1:N
            xBestValue(i) = func(xBest(i,:));
        end
        temp = find(xBestValue == min(xBestValue));
        gBest = xBest(temp(1),:);              %Ⱥ������λ��
        gBestValue = xBestValue(temp(1));      %Ⱥ������λ�õ���Ӧ��
        xBest(temp(1),:) = xBest(1,:);
        xBestValue(temp(1)) = xBestValue(1);
        xBest(1,:) = gBest;
        xBestValue(1) = gBestValue;
        
        %%%%%%%%%%%%%%%%%%%% ʨȺ��ɫ���� %%%%%%%%%%%%%%%%%%%%
        %ʨ��λ�ó�ʼ��
        xKing = xBest(1,:);
        xKingValue = xBestValue(1);
        %ĸʨ��λ�ó�ʼ��
        n = floor(beta * N);                   %����ʨ������
        xHunter = xBest(2:n,:);
        xHunterValue = xBestValue(2:n);
        %ĸʨ����λ�ó�ʼ��
        temp = find(xHunterValue == min(xHunterValue));
        xHunterBest = xHunter(temp(1),:);
        %��ʨλ�ó�ʼ��
        xCub = xBest(n+1:N,:);
        xCubValue = xBestValue(n+1:N);
        gWorst = zeros(1,D);
        
        %%%%%%%%%%%%%%%%%%%% �������� %%%%%%%%%%%%%%%%%%%%
        for i = 1:T
            
            %%%%%%%%%%%% ����һ�������¼ʨȺ·�� %%%%%%%%%%%%%%
            if k == C
                %��ɢ��
                hunterDens = mean(var(xHunter)) / ((xMax-xMin)^2/12);
                densTrac(i) = hunterDens;
                %��¼ȫ�����ŵ�仯·��
                xKingTrace(i,:) = gBest;     
                
                if trac == 1 && D == 2 && rem(i,5)-1 == 0 && i < 30
                    subplot(2,3,s);
                    %���Ƶȸ���
                    drawcontour(func,xMax,xMin);
                    hold on
                    %���Ƹ�����ʷ���ŷֲ�
                    scatter(xHunter(1:n-1,1), xHunter(1:n-1,2),10,'filled','bo');
                    hold on
                    scatter(xCub(1:N-n,1), xCub(1:N-n,2),10,'filled','ko');
                    hold on        
                    scatter(xBest(1,1), xBest(1,2),10,'filled','ro');
                    title(['t=',num2str(i),' ��ɢ��=',num2str(roundn(hunterDens,-2))]);     
                    grid on
                    hold off
                    s = s + 1;
                elseif trac == 2 && D == 2
                    %���Ƶȸ���
                    drawcontour(func,xMax,xMin);
                    hold on
                    %��ǰ��Ⱥ�ֲ�                  
                    scatter(xHunter(1:n-1,1), xHunter(1:n-1,2),10,'filled','bo');
                    hold on
                    scatter(xCub(1:N-n,1), xCub(1:N-n,2),10,'filled','ko');
                    hold on        
                    scatter(xKing(1,1), xKing(1,2),10,'filled','ro');
                    title(['t=',num2str(i),' ��ɢ��=',num2str(roundn(hunterDens,-2))]);     
                    grid on
                    hold off

                    frame = getframe;
                    frame.cdata = imresize(frame.cdata, [435 343]); %������Ƶ��ߣ�HΪ����(��)��WΪ����(��)
                    writeVideo(aviObj,frame);
                end  
            end
            
            %%%%%%%%%%%%%%%%%%%% ����λ�� %%%%%%%%%%%%%%%%%%%%
            alphaF = step * exp(-30*((i/T)^10));    %ĸʨ�ƶ���Χ�Ŷ�����
            alphaC = step * (T - i) / T;            %��ʨ�ƶ���Χ�Ŷ�����
            %����ʨ��λ��
            xKing = gBest * (1 + rand * norm(xKing - gBest));
            %����ĸʨλ��
            for j = 1:n-1
                t = ceil(rand*(n-1));                   %���ѡ��Э�����
                while t == j
                    t = ceil(rand*(n-1));
                end   
                xHunter(j,:) = ((xBest(j+1,:)+xBest(t+1,:))/2) * (1 + alphaF * rand);
                %�߽���������
                for t = 1:D
                    if (xHunter(j,t) > xMax) || (xHunter(j,t) < xMin)
                        xHunter(j,t) = rand * (xMax - xMin) + xMin;
                    end
                end
            end
            %������ʨλ��	
%             gWorst = mean(min(xCub)) + mean(max(xCub)) - gBest;
            gWorst = xMax + xMin -gBest;
            
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

            %%%%%%%%%%%%%%%%%%%% ������ʷ���ż�¼ %%%%%%%%%%%%%%%%%%%%
            %����ʨ����ʷ����
            xKingValue = func(xKing);
            if xKingValue < xBestValue(1)
                xBest(1,:) = xKing;
                xBestValue(1) = xKingValue;
            end
            %����ĸʨ��ʷ����
            for t = 1:n-1
                xHunterValue(t) = func(xHunter(t,:));   
                if xHunterValue(t) < xBestValue(t+1)
                    xBest(t+1,:) = xHunter(t,:);
                    xBestValue(t+1) = xHunterValue(t);
                end
            end
            %������ʨ��ʷ����
            for t = 1:N-n
                xCubValue(t) = func(xCub(t,:));
                if xCubValue(t) < xBestValue(n + t)
                    xBest(n + t,:) = xCub(t,:);
                    xBestValue(n + t) = xCubValue(t);
                end
            end
            %����xHunterBest
            temp = find(xHunterValue == min(xHunterValue));
            xHunterBest = xHunter(temp(1),:);
            
            %����gBest��ÿ����10������ȷ��ʨ��λ��
            temp = find(xBestValue == min(xBestValue));
            gBest = xBest(temp(1),:);              
            gBestValue = xBestValue(temp(1));      
            %gBestÿ�ε������п��ܸ��£����xKingֻ�ڳ�ʼ�׶�һ������gBest
            if rem(i,10) == 0 && gBestValue ~= xBestValue(1)
                xBest(temp(1),:) = xBest(1,:);
                xBestValue(temp(1)) = xBestValue(1);
                xBest(1,:) = gBest;
                xBestValue(1) = gBestValue;
            end    
            
            gbtemp(k,i) = gBestValue;
        end
    end

    gb = mean(gbtemp);
    
    %%%%%%%%%%%%%%%%%%%% ������ %%%%%%%%%%%%%%%%%%%%
    if trac == 2 && D == 2
        close(aviObj);
    end
    
    %��Ⱥ�ֲ��ط���
    if trac ~= 0  
        figure
        plot(densTrac) 
        title('ĸʨ��ɢ��')
    end
    
    %ȫ�����ŵ�Ľ���·��
    if trac == 1 && D == 2
        figure
        %���Ƶȸ���
        drawcontour(func,xMax,xMin);
        hold on
        %��������ʨȺ�ֲ������һ�����飩
        scatter(xBest(2:n,1), xBest(2:n,2),10,'filled','bo');
        hold on
        scatter(xBest(n+1:N,1), xBest(n+1:N,2),10,'filled','ko');
        hold on
        %ȫ�����ŵ�仯
        plot(xKingTrace(:,1), xKingTrace(:,2),'ro-');
        title('ȫ�����ŵ����·��')
        grid on
        hold off
    end

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

%���Ƶȸ���
function [] = drawcontour(func,xMax,xMin)
    x = linspace(xMin,xMax,100);
    y = linspace(xMin,xMax,100);
    zz = zeros(100,100);
    [xx,yy] = meshgrid(x,y);
    for t = 1:100
        for j = 1:100
            zz(t,j) = func([xx(t,j),yy(t,j)]);
        end
    end
    contour(xx,yy,zz,10);
end