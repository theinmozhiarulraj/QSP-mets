npatients = length(params_out.patient);
ntumors=length(TumorArray);
x_time=simDataPSA(params_out.iPatient(1)).simData.Time;
tumornames={'Primary tumor','Lung Met','','Other Met'};
C={'y',[0.4660 0.6740 0.1880],[0.4660 0.6740 0.1880],[0 0.4470 0.7410]};


figure();

for k=1:4

    if(k==1 || k==2 ||k==4)

        disp(k)
        count=1;
        % change in tumor dia at the end of treatment
        tumor_timeseries=NaN(length(x_time),1);

    end

    for i = 1:npatients

       if(params_out.patient(i)==1)
                
                pos = find(ismember(simDataPSApost(i).simData.DataNames,['D_T' TumorArray{k}]));
                D_T_temp = simDataPSApost(i).simData.Data(:,pos);
    
                if(D_T_temp(1)>=0.2)
    
                    pos_perc = find(ismember(simDataPSApost(i).simData.DataNames,['D_T' TumorArray{k} '_perc']));
                    D_T_temp_perc = simDataPSApost(i).simData.Data(:,pos_perc);
    
                    tumor_timeseries(:,count)=D_T_temp_perc;

                    count=count+1;
    
                end
    
        end

    end

    if(k==1 || k==2 ||k==4)

        % row wise median
        med_tumorseries=median(tumor_timeseries,2);
    
        % confidence intervals
        N_tumorseries=size(tumor_timeseries,1);
        ySEM_tumorseries=std(tumor_timeseries,[],2)/sqrt(N_tumorseries);
        CI95_tumorseries=tinv([0.025 0.975],N_tumorseries-1);
        yCI95_tumorseries=med_tumorseries+CI95_tumorseries.*ySEM_tumorseries;
    
        % plot here
        plot(x_time,med_tumorseries,'-','color',C{k},'DisplayName',tumornames{k},'LineWidth',3);
        patch([x_time' fliplr(x_time')], [yCI95_tumorseries(:,1)' fliplr(yCI95_tumorseries(:,2)')], C{k}, 'FaceAlpha',0.1, 'EdgeColor','none')
        hold on;

    end

end

xlabel('Time (days)', 'FontSize', 20)
ylabel('% change in tumor size', 'FontSize', 20)

ax=gca;
ax.FontSize = 20;
xlim([0 400])
%title(name, 'FontSize', 16)
%legend([],'FontSize',18)



