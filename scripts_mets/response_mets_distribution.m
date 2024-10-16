% organs=["primary","lung","others"];
% suffix=["","_Ln1","_other"];

npatients = length(params_out.patient);
ntumors=length(TumorArray);

% change in tumor dia at the end of treatment
primary=NaN(npatients*2,1);
lung=NaN(npatients*2,1);
others=NaN(npatients*2,1);

% change in tumor dia at best response
primary_best=NaN(npatients*2,1);
lung_best=NaN(npatients*2,1);
others_best=NaN(npatients*2,1);

nprimary=0;
nlung=0;
nothers=0;


for i = 1:npatients

   if(params_out.patient(i)==1)

       for j=1:ntumors
            
            pos = find(ismember(simDataPSApost(i).simData.DataNames,['D_T' TumorArray{j}]));
            D_T_temp = simDataPSApost(i).simData.Data(:,pos);

            if(D_T_temp(1)>=0.2)

                pos_perc = find(ismember(simDataPSApost(i).simData.DataNames,['D_T' TumorArray{j} '_perc']));
                D_T_temp_perc = simDataPSApost(i).simData.Data(:,pos_perc);

                if(TumorArray{j}=="")

                    primary(i,1)=D_T_temp_perc(end);
                    primary_best(i,1)=min(D_T_temp_perc);
                    nprimary=nprimary+1;

                elseif(TumorArray{j}=="_Ln1")

                    lung((i*2)-1,1)=D_T_temp_perc(end);
                    lung_best((i*2)-1,1)=min(D_T_temp_perc);
                    nlung=nlung+1;

                elseif(TumorArray{j}=="_Ln2")

                    lung(i*2,1)=D_T_temp_perc(end);                    
                    lung_best(i*2,1)=min(D_T_temp_perc);
                    nlung=nlung+1;


                elseif(TumorArray{j}=="_other")

                    others(i,1)=D_T_temp_perc(end);
                    others_best(i,1)=min(D_T_temp_perc);         
                    nothers=nothers+1;

                else

                    %dont do anything

                end

            end

       end
       
   end

end

figure();
histogram(primary,10,'FaceAlpha',0.2,'Normalization','pdf')
hold on;
histogram(lung,10,'FaceAlpha',0.6,'Normalization','pdf')
hold on;
histogram(others,10,'FaceAlpha',0.2,'Normalization','pdf')

dataset=table(primary,lung,others);
writetable(dataset,'response_mets_end.csv')

dataset_best=table(primary_best,lung_best,others_best);
writetable(dataset_best,'response_mets_best.csv')