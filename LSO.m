function [sucRate,minv,maxv,meanv,stdv,convTime] = LSO(func,xMax,xMin,target,error,C,T,D,N,beta,trac,fitness)
    %输入参数
    %func:测试函数
    %xMax:函数定义域
    %target:寻优目标
    %error:允许误差
    %C:试验次数
    %T:最大迭代次数
    %D:函数维度
    %N:种群数量
    %beta:成年狮所占比例因子
    %trac:是否绘制二维种群路径
    %fitness:是否绘制适应度进化曲线
    
    %返回值
    %sucRate:成功率
    %min:寻优结果最小值
    %max:寻优结果最大值
    %mean:寻优结果平均值
    %std:寻优结果标准差
    %convTime:收敛所需迭代次数
    
    step = 0.05 * (xMax - xMin);

    gbtemp = zeros(C,T);                %每次测试、每次迭代的结果
    xKingTrace = zeros(T,D);
    densTrac =zeros(1,T);

    if trac == 1 && D == 2
        figure
        s = 1;  %子图
    elseif trac == 2 && D == 2
        figure
        aviObj = VideoWriter('AdaLSO.avi');%保存为avi
        aviObj.FrameRate = 20;
        open(aviObj);
    end
    %%%%%%%%%%%%%%%%%%%% 开始试验 %%%%%%%%%%%%%%%%%%%%
    for k = 1:C

        %%%%%%%%%%%%%%%%%%%% 狮群状态初始化 %%%%%%%%%%%%%%%%%%%%
        %初始化 个体历史最优位置xBest 和 群体历史最优位置gBest
        xInicial = rand(N,D) * (xMax - xMin) + xMin;
        xBest = xInicial;                   
        xBestValue = ones(1,N);             %个体最优位置的适应度
        for i = 1:N
            xBestValue(i) = func(xBest(i,:));
        end
        temp = find(xBestValue == min(xBestValue));
        gBest = xBest(temp(1),:);              %群体最优位置
        gBestValue = xBestValue(temp(1));      %群体最优位置的适应度
        xBest(temp(1),:) = xBest(1,:);
        xBestValue(temp(1)) = xBestValue(1);
        xBest(1,:) = gBest;
        xBestValue(1) = gBestValue;
        
        %%%%%%%%%%%%%%%%%%%% 狮群角色分配 %%%%%%%%%%%%%%%%%%%%
        %狮王位置初始化
        xKing = xBest(1,:);
        xKingValue = xBestValue(1);
        %母狮的位置初始化
        n = floor(beta * N);                   %成年狮的数量
        xHunter = xBest(2:n,:);
        xHunterValue = xBestValue(2:n);
        %母狮最优位置初始化
        temp = find(xHunterValue == min(xHunterValue));
        xHunterBest = xHunter(temp(1),:);
        %幼狮位置初始化
        xCub = xBest(n+1:N,:);
        xCubValue = xBestValue(n+1:N);
        gWorst = zeros(1,D);
        
        %%%%%%%%%%%%%%%%%%%% 迭代过程 %%%%%%%%%%%%%%%%%%%%
        for i = 1:T
            
            %%%%%%%%%%%% 抽样一次试验记录狮群路径 %%%%%%%%%%%%%%
            if k == C
                %离散度
                hunterDens = mean(var(xHunter)) / ((xMax-xMin)^2/12);
                densTrac(i) = hunterDens;
                %记录全局最优点变化路径
                xKingTrace(i,:) = gBest;     
                
                if trac == 1 && D == 2 && rem(i,5)-1 == 0 && i < 30
                    subplot(2,3,s);
                    %绘制等高线
                    drawcontour(func,xMax,xMin);
                    hold on
                    %绘制个体历史最优分布
                    scatter(xHunter(1:n-1,1), xHunter(1:n-1,2),10,'filled','bo');
                    hold on
                    scatter(xCub(1:N-n,1), xCub(1:N-n,2),10,'filled','ko');
                    hold on        
                    scatter(xBest(1,1), xBest(1,2),10,'filled','ro');
                    title(['t=',num2str(i),' 分散度=',num2str(roundn(hunterDens,-2))]);     
                    grid on
                    hold off
                    s = s + 1;
                elseif trac == 2 && D == 2
                    %绘制等高线
                    drawcontour(func,xMax,xMin);
                    hold on
                    %当前种群分布                  
                    scatter(xHunter(1:n-1,1), xHunter(1:n-1,2),10,'filled','bo');
                    hold on
                    scatter(xCub(1:N-n,1), xCub(1:N-n,2),10,'filled','ko');
                    hold on        
                    scatter(xKing(1,1), xKing(1,2),10,'filled','ro');
                    title(['t=',num2str(i),' 分散度=',num2str(roundn(hunterDens,-2))]);     
                    grid on
                    hold off

                    frame = getframe;
                    frame.cdata = imresize(frame.cdata, [435 343]); %设置视频宽高：H为行数(高)，W为列数(宽)
                    writeVideo(aviObj,frame);
                end  
            end
            
            %%%%%%%%%%%%%%%%%%%% 更新位置 %%%%%%%%%%%%%%%%%%%%
            alphaF = step * exp(-30*((i/T)^10));    %母狮移动范围扰动因子
            alphaC = step * (T - i) / T;            %幼狮移动范围扰动因子
            %更新狮王位置
            xKing = gBest * (1 + rand * norm(xKing - gBest));
            %更新母狮位置
            for j = 1:n-1
                t = ceil(rand*(n-1));                   %随机选择协作伙伴
                while t == j
                    t = ceil(rand*(n-1));
                end   
                xHunter(j,:) = ((xBest(j+1,:)+xBest(t+1,:))/2) * (1 + alphaF * rand);
                %边界条件处理
                for t = 1:D
                    if (xHunter(j,t) > xMax) || (xHunter(j,t) < xMin)
                        xHunter(j,t) = rand * (xMax - xMin) + xMin;
                    end
                end
            end
            %更新幼狮位置	
%             gWorst = mean(min(xCub)) + mean(max(xCub)) - gBest;
            gWorst = xMax + xMin -gBest;
            
            for j = 1:N-n
                q = rand;
                if q <= 1/3     %向狮王靠近
                    xCub(j,:) = ((gBest + xBest(j + n,:))/2) * (1 + alphaC * rand);
                elseif q <2/3	%向母狮中的最优位置xHunterBest靠近
                    xCub(j,:) = ((xHunterBest + xBest(j + n,:))/2) * (1 + alphaC * rand);
                else            %远离狮王
                    xCub(j,:) = ((gWorst + xBest(j + n,:))/2) * (1 + alphaC * rand);
                end
                %边界条件处理
                for t = 1:D
                    if (xCub(j,t) > xMax) || (xCub(j,t) < xMin)
                        xCub(j,t) = rand * (xMax - xMin) + xMin;
                    end
                end
            end

            %%%%%%%%%%%%%%%%%%%% 更新历史最优记录 %%%%%%%%%%%%%%%%%%%%
            %更新狮王历史最优
            xKingValue = func(xKing);
            if xKingValue < xBestValue(1)
                xBest(1,:) = xKing;
                xBestValue(1) = xKingValue;
            end
            %更新母狮历史最优
            for t = 1:n-1
                xHunterValue(t) = func(xHunter(t,:));   
                if xHunterValue(t) < xBestValue(t+1)
                    xBest(t+1,:) = xHunter(t,:);
                    xBestValue(t+1) = xHunterValue(t);
                end
            end
            %更新幼狮历史最优
            for t = 1:N-n
                xCubValue(t) = func(xCub(t,:));
                if xCubValue(t) < xBestValue(n + t)
                    xBest(n + t,:) = xCub(t,:);
                    xBestValue(n + t) = xCubValue(t);
                end
            end
            %更新xHunterBest
            temp = find(xHunterValue == min(xHunterValue));
            xHunterBest = xHunter(temp(1),:);
            
            %更新gBest，每迭代10次重新确定狮王位置
            temp = find(xBestValue == min(xBestValue));
            gBest = xBest(temp(1),:);              
            gBestValue = xBestValue(temp(1));      
            %gBest每次迭代都有可能更新，因此xKing只在初始阶段一定等于gBest
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
    
    %%%%%%%%%%%%%%%%%%%% 输出结果 %%%%%%%%%%%%%%%%%%%%
    if trac == 2 && D == 2
        close(aviObj);
    end
    
    %种群分布熵分析
    if trac ~= 0  
        figure
        plot(densTrac) 
        title('母狮分散度')
    end
    
    %全局最优点的进化路径
    if trac == 1 && D == 2
        figure
        %绘制等高线
        drawcontour(func,xMax,xMin);
        hold on
        %绘制最终狮群分布（最后一次试验）
        scatter(xBest(2:n,1), xBest(2:n,2),10,'filled','bo');
        hold on
        scatter(xBest(n+1:N,1), xBest(n+1:N,2),10,'filled','ko');
        hold on
        %全局最优点变化
        plot(xKingTrace(:,1), xKingTrace(:,2),'ro-');
        title('全局最优点更新路径')
        grid on
        hold off
    end

    %数据分析：成功率、误差范围
    sucRate = 0;
    convTime = 0;
    for i = 1:C
        if gbtemp(i,T) - target < error     %达到迭代次数上限时满足要求则成功
            sucRate = sucRate + 1;
        end
    end
    sucRate = sucRate / C;
    minv = min(gbtemp(:,T));
    maxv = max(gbtemp(:,T));
    meanv = mean(gbtemp(:,T));
    stdv = std(gbtemp(:,T));
    
    %收敛所需（平均）迭代次数
    for i = 1:T
        if gb(i) - target > error
            convTime = convTime + 1;
        end
    end  
    
    %适应度进化曲线
    if fitness == 1
        figure
        plot(gb)
        xlabel('迭代次数');
        ylabel('适应度值f(x)');
        title(['适应度进化曲线,维度D=',num2str(D)])
    end
    
end

%绘制等高线
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