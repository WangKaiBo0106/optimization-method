function [sucRate,minv,maxv,meanv,stdv,convTime] = PSO(func,xMax,xMin,target,error,C,T,D,N,~,~,fitness)
    c1 = 2;
    c2 = 2;
    Wmax = 0.8;
    Wmin = 0.4;
%     w = 0.7298;

    Vmax = 20;
    Vmin = -10;
    gbtemp = zeros(C,T);                %每次测试、每次迭代的结果
    for k = 1:C
        %初始化种群个体
        x = rand(N,D) * (xMax - xMin) + xMin;	%随机生成N个*D维
        v = rand(N,D) * (Vmax - Vmin) + Vmin;

        %初始化个体最优位置和最优值
        p = x;			%p--个体最优位置，N个
        pbest = ones(N,1);	%pbest--个体最优值，N个
        for i = 1:N
            pbest(i) = func(x(i,:));
        end
        %初始化全局最优位置和最优值
        g = ones(1,D);		%g--全局最优位置，1个
        gbest = inf;		%gbest--全局最优值，1个
        for i = 1:N
            if(pbest(i) < gbest)
                g = p(i,:);
            gbest = pbest(i);
            end
        end

        %按照公式迭代
        for i = 1:T
            for j = 1:N
            %更新个体最优位置和最优值
            if (func(x(j,:)) < pbest(j))
                p(j,:) = x(j,:);
                pbest(j) = func(x(j,:));
            end
            %更新全局最优位置和最优值
            if (pbest(j) < gbest)
                g = p(j,:);
                gbest = pbest(j,:);
            end

            %计算动态惯性权重
            w = Wmax - (Wmax -Wmin)*i/T; %逐渐减小

            %更新位置和速度值
            v(j,:) = w*v(j,:) + c1*rand*(p(j,:)-x(j,:)) + c2*rand*(g-x(j,:));
            x(j,:) = x(j,:) + v(j,:);
            %边界条件处理
            for ii = 1:D
                if (v(j,ii) > Vmax) || (v(j,ii) < Vmin)
                v(j,ii) = rand * (Vmax - Vmin) + Vmin;
                end
                if (x(j,ii) > xMax) || (x(j,ii) < xMin)
                x(j,ii) = rand * (xMax - xMin) + xMin;
                end
            end
            end
            %记录历代全局最优值
            gbtemp(k,i) = gbest;
        end
    end
    
    gb = mean(gbtemp);
    
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


