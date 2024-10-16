%thein: this script plots the response of individual tumors from randomly
%sampled patients

% response status
response_categories={'CR/PR','SD'};

% colors
C={'k','b',[0.3010, 0.7450, 0.9330],[0.3010, 0.7450, 0.9330],[0.4660, 0.6740, 0.1880]};

% index of 1st patient
index_1patient=params_out.iPatient(1);

% tumor names
display_names={'All tumors','Primary','Lung met','Lung met','Other met'};
names_percentchange={'_all',TumorArray{1:end}};
tumorperc_pos=zeros(length(names_percentchange),1);
tumordia_pos=zeros(length(names_percentchange),1);

for k=1:length(names_percentchange)
    tumorperc_pos(k,1) = [find(strcmp(simDataPSApost(index_1patient).simData.DataNames,['D_T' names_percentchange{k} '_perc']) )];
    tumordia_pos(k,1) = [find(strcmp(simDataPSApost(index_1patient).simData.DataNames,['D_T' names_percentchange{k}]) )];
end

%time
time_temp=simDataPSA(index_1patient).simData.Time;

% number of patients to sample for each response category 
nsamples=10;

for i=1:length(response_categories)
    % get the index of all patients with given response status
    index_patients=find(ismember(params_out.RECIST_mets, response_categories{i})==1);

    % randomly sample patients
    index_sampled=randsample(index_patients,nsamples);

    for j=1:nsamples

        figure();
        temp_all_tumors=simDataPSApost(index_sampled(j)).simData.Data(:,tumorperc_pos(1,1));
        plot(time_temp,temp_all_tumors,'-','Color','k','LineWidth',2.5,'DisplayName','All tumors');
        hold on;

        % plot all tumors change in diameter
        for k=2:length(names_percentchange)
            
            temp_tumordia=simDataPSApost(index_sampled(j)).simData.Data(:,tumordia_pos(k,1));

            if(temp_tumordia(1)>=0.2)
                temp=simDataPSApost(index_sampled(j)).simData.Data(:,tumorperc_pos(k,1));
                plot(time_temp,temp,':','Color',C{k},'LineWidth',2.5,'DisplayName',display_names{k});
            end
            hold on;
        
        end

        ylim([-100 100])
        ylabel('Percent change in tumor diameter','Fontsize',14)
        xlabel('Time (days)','Fontsize',14)
        legend();
        title(response_categories{i},'Fontsize',14)
        set(gca,'Fontsize',14)
        


    end

end

