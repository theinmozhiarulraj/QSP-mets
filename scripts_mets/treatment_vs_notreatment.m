% load the data from untreated and treated 
simdata_untreated = load('simDataPSA_v5_reference.mat').simDataPSA;
simdata_treated = load('simDataPSA_v5.mat').simDataPSA;
params_out_untreated=load('params_out_v5.mat').params_out;

% note: params_out should also be loaded for this script to work
patientno=params_out_untreated.iPatientPlaus(1);
npatients=length(params_out_untreated.RECIST);

% number of species
[~,nspecies]=size(simdata_untreated(patientno).simData.Data);
[~,nspecies_test]=size(simdata_treated(patientno).simData.Data);
ntimepts=length(simdata_untreated(patientno).simData.Time);
timepoints=simdata_untreated(patientno).simData.Time;

if(nspecies~=nspecies_test)
    disp("ERROR: number of species in both datasets do not match")
end

nper_col=3;
nper_row=3;
nper_fig=nper_col*nper_row;

start=1;
for j=1:(nspecies/nper_fig)+1

    figure()
    tiledlayout(nper_row,nper_col);

    for i=start:start+nper_fig-1

        difference_crpr=NaN(npatients,ntimepts);
        difference_sd=NaN(npatients,ntimepts);
        difference_pd=NaN(npatients,ntimepts);

        start=start+nper_fig; 

        compartment=simdata_untreated(1).simData.DataInfo{i,1}.Compartment;
        species=simdata_treated(1).simData.DataInfo{i,1}.Name;

        for k=1:npatients

            if(params_out_untreated.patient(k)==1)
        
                if(params_out_untreated.RECIST_mets(k)=="CR/PR")
        
                   difference_crpr(k,:)=(simdata_treated(k).simData.Data(:,i)-simdata_untreated(k).simData.Data(:,i))*100./simdata_untreated(k).simData.Data(:,i);
        
                elseif(params_out_untreated.RECIST_mets(k)=="SD")
        
                   difference_sd(k,:)=(simdata_treated(k).simData.Data(:,i)-simdata_untreated(k).simData.Data(:,i))*100./simdata_untreated(k).simData.Data(:,i);
        
                elseif(params_out_untreated.RECIST_mets(k)=="PD")
        
                   difference_pd(k,:)=(simdata_treated(k).simData.Data(:,i)-simdata_untreated(k).simData.Data(:,i))*100./simdata_untreated(k).simData.Data(:,i);
        
                end

            end

        end

        med_crpr=median(difference_crpr,'omitnan');
        med_sd=median(difference_sd,'omitnan');
        med_pd=median(difference_pd,'omitnan');
    
        nexttile
        plot(timepoints,med_crpr,'DisplayName','CR/PR','LineWidth',2);
        hold on;
        plot(timepoints,med_sd,'DisplayName','SD','LineWidth',2);
        hold on;
        plot(timepoints,med_pd,'DisplayName','PD','LineWidth',2);        
        
        xlabel("Time (days)", 'FontSize', 14);
        ylabel("Difference (Treated-untreated)", 'FontSize', 14);

        ax=gca;
        ax.FontSize = 13;
    
        xlim([0 100])
        title([species ' ' compartment]);
        legend('FontSize',13)

    end

end


