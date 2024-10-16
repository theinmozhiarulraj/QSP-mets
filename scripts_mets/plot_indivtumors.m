% plot response of different tumors side by side

npatients = length(params_out.iPatient);
ntumors=length(TumorArray);

%tumordia=NaN(npatients,ntumors);
tumordia=NaN(20,ntumors);

%for i = 1:npatients
for i = 1:20

    for j=1:ntumors
        index = params_out.iPatient(i);
        pos = find(ismember(simDataPSApost(index).simData.DataNames,['D_T' TumorArray{j} '_perc']));
        D_T_temp = simDataPSApost(index).simData.Data(:,pos);
        tumordia(i,j)=D_T_temp(length(D_T_temp));
    end

end
figure();
bar(tumordia)

[rows columns]=size(tumordia);
for j=1:rows
    for k=1:columns

      if(tumordia(j,k)==0)
        tumordia(j,k)=NaN;
      end

        if(tumordia(j,k)==-0.0000)
         tumordia(j,k)=NaN;
        end

    end
end
figure();
scatter(1:20,tumordia,'filled')



% sample n patients with different response status and plot them
nsamples=10;

% find the index of patients with different response status
patients_CRPR=find(params_out.RECIST_mets=="CR/PR");
patients_SD=find(params_out.RECIST_mets=="SD");
patients_PD=find(params_out.RECIST_mets=="PD");

% randomly permut the indices
indices = randperm(length(patients_CRPR));  
patients_CRPR_final = patients_CRPR(indices);

indices = randperm(length(patients_SD));  
patients_SD_final = patients_SD(indices);

indices = randperm(length(patients_PD));  
patients_PD_final = patients_PD(indices);

nsamples_tmp_CRPR=nsamples;
nsamples_tmp_SD=nsamples;
nsamples_tmp_PD=nsamples;

if(length(patients_CRPR_final)<nsamples)
    nsamples_tmp_CRPR=length(patients_CRPR_final);
end

if(length(patients_SD_final)<nsamples)
    nsamples_tmp_SD=length(patients_SD_final);
end

if(length(patients_PD_final)<nsamples)
    nsamples_tmp_PD=length(patients_PD_final);
end

ntumors=length(TumorArray);

tumordia_CRPR=NaN(nsamples_tmp_CRPR,ntumors);
tumordia_SD=NaN(nsamples_tmp_SD,ntumors);
tumordia_PD=NaN(nsamples_tmp_PD,ntumors);


for j=1:ntumors

        for i=1:nsamples_tmp_CRPR
            pos = find(ismember(simDataPSApost(patients_CRPR_final(i)).simData.DataNames,['D_T' TumorArray{j} '_perc']));
            D_T_temp = simDataPSApost(patients_CRPR_final(i)).simData.Data(:,pos);
            tumordia_CRPR(i,j)=D_T_temp(length(D_T_temp));
        end

        for i=1:nsamples_tmp_SD
            pos = find(ismember(simDataPSApost(patients_SD_final(i)).simData.DataNames,['D_T' TumorArray{j} '_perc']));
            D_T_temp = simDataPSApost(patients_SD_final(i)).simData.Data(:,pos);
            tumordia_SD(i,j)=D_T_temp(length(D_T_temp));
        end

         for i=1:nsamples_tmp_PD
            pos = find(ismember(simDataPSApost(patients_PD_final(i)).simData.DataNames,['D_T' TumorArray{j} '_perc']));
            D_T_temp = simDataPSApost(patients_PD_final(i)).simData.Data(:,pos);
            tumordia_PD(i,j)=D_T_temp(length(D_T_temp));
        end


end

figure();
bar(tumordia_CRPR)

figure();
bar(tumordia_SD)

figure();
bar(tumordia_PD)

