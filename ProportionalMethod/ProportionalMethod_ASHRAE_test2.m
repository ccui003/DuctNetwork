% Starting from the second lowest P
ASHRAE_duct;
DesignFlow = {[700,250,950],[275,275,475,475,200,200]};% M
lambda = 0.5;

N = cellfun(@length,DesignFlow);
Branch = length(N);% M, number of branches
Q = cell(1,Branch);
P = cell(1,Branch);
ERR = cell(1,Branch);
I = cell(1,Branch);
TargetFlow = cell(1,Branch);
Theta = cell(1,Branch);
duct.S(206) = 500; % Fan Max Pressure;
duct.S(207) = 3000; % Fan Max Flow;
duct.SetDamperAndReadFlow(zeros(1,sum(N)),1:sum(N));% fully open all dampers
options = optimset('Display','none','TolX',0.1);

for i = 1:Branch
    DamperID = sum(N(1:i-1))+1 : sum(N(1:i));% dampers being adjusted in the current branch
    Q{i} = duct.SetDamperAndReadFlow([],DamperID);
    P{i} = Q{i}./DesignFlow{i};
    ERR{i} = abs(P{i}-1);
    Iter = 0;
    Theta{i} = zeros(0,N(i));
    while max(ERR{i}(end,:))>0.1 && Iter <50
        [~,I{i}(end+1,:)] = sort(P{i}(end,:));
        TargetFlow{i}(end+1,:) = lambda*DesignFlow{i}*P{i}(end,I{i}(end,1))+(1-lambda)*Q{i}(end,:);
        Iter = size(Theta{i},1)+1;
        if Iter>1
            Theta{i}(Iter,I{i}(end,1)) = Theta{i}(Iter-1,I{i}(end,1));
        end
        fprintf('Start to adjust dampers, Iter = %d\n',Iter);
        for j = I{i}(end,2:end)
            theta = fminbnd(@(theta)abs(duct.SetDamperAndReadFlow(theta,DamperID(j))-TargetFlow{i}(end,j)),0,90,options);
            Theta{i}(Iter,j) = theta;
            duct.SetDamperAndReadFlow(theta,DamperID(j)); 
        end
        theta = Theta{i}(end,:)
        Theta{i}(end+1,:)=Theta{i}(end,:)-min(Theta{i}(end,:));
        theta = Theta{i}(end,:)
        Q{i}(end+1,:) = duct.SetDamperAndReadFlow(Theta{i}(end,:),DamperID);
        q = Q{i}(end,:)
        P{i}(end+1,:) = Q{i}(end,:)./DesignFlow{i};
        p = P{i}(end,:)
        ERR{i}(end+1,:) = abs(P{i}(end,:)-1);
        error = ERR{i}(end,:)
        q_fan = sum(Q{i}(end,:));
        T_fan = sum(DesignFlow{i});
        if q_fan>1.05*T_fan
            r = T_fan/q_fan;
        elseif q_fan<0.95*T_fan
            r = T_fan/q_fan;
        else
            r=1;
        end
        r
        duct.S(206) = duct.S(206)*r^2; % Fan Max Pressure;
%         Fan_P = str2double(get_param([SystemName,'/Fan'],'Max_P'));
%         set_param([SystemName,'/Fan'],'Max_P',num2str(Fan_P*r^2));
        duct.S(207) = duct.S(207)*r; % Fan Max Flow;
%         Fan_q = str2double(get_param([SystemName,'/Fan'],'Max_q'));
%         set_param([SystemName,'/Fan'],'Max_q',num2str(Fan_q*r)); 
      
        Q{i}(end+1,:) = duct.SetDamperAndReadFlow([],DamperID);
        q = Q{i}(end,:)
        P{i}(end+1,:) = Q{i}(end,:)./DesignFlow{i};
        p = P{i}(end,:)
        ERR{i}(end+1,:) = abs(P{i}(end,:)-1);
        error = ERR{i}(end,:)
    end
end

Q_final = duct.SetDamperAndReadFlow([],1:sum(N));
P_final = Q_final./[DesignFlow{:}];
ERR_final = abs(P_final-1);
FileName = ['ASHRAE','_ProportionalMethod_',datestr(now,'ddmmyyyyHHMMSS')];
save([FileName,'.mat'],'Q','P','ERR','I','TargetFlow','Theta','Q_final','P_final','ERR_final');

%%
for i = 1:Branch
    figure(i)
    set(gcf,'Position',[130 119 1538 858])
    clf
    subplot(2,2,1)
    plot(Q{i});
    hold on
    plot(ones(size(Q{i},1),1)*DesignFlow{i},'-k')
    Tag = arrayfun(@(damper)sprintf('Damper %d',damper),sum(N(1:i-1))+1 : sum(N(1:i)),'UniformOutput',false);
    Tag{end+1}='Design Flow';
    legend(Tag);
    ylabel('Flow Rate (l/s)')
    subplot(2,2,2)
    plot(P{i}*100)
    hold on
    plot(ones(size(Q{i},1),1)*100,'-k')
    plot(ones(size(Q{i},1),1)*[0.9,1.1]*100,'--k')
    ylabel('Ratio to Design Flow Rate (%)')
    subplot(2,2,3)
    plot(TargetFlow{i})
    ylabel('Target Flow Rate (l/s)')
    subplot(2,2,4)
    plot(Theta{i})
    ylabel('Damper Angle (degree)')
    savefig([FileName,'_Branch',num2str(i),'.fig'])
    saveas(gcf,[FileName,'_Branch',num2str(i),'.png'])
end
