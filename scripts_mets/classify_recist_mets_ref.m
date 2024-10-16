% calculate percent patients with PD, SD, CR, PR and PR/CR
% only plausible patients

npatients_plausible = length(params_out.iPlaus);

calc_PD=0;
calc_CR=0;
calc_PR=0;
calc_SD=0;
calc_NP=0; % not patient
vector_cr=zeros(1,npatients_plausible);
vector_pr=zeros(1,npatients_plausible);
vector_sd=zeros(1,npatients_plausible);
vector_pd=zeros(1,npatients_plausible);

for i=1:npatients_plausible

    j=params_out.iPlaus(i);

    if(params_out.RECIST_mets_ref{j}=="CR")
        calc_CR=calc_CR+1;
        vector_cr(i)=1;
    elseif(params_out.RECIST_mets_ref{j}=="PR")
        calc_PR=calc_PR+1;
        vector_pr(i)=1;        
    elseif(params_out.RECIST_mets{j}=="SD")
        calc_SD=calc_SD+1;
        vector_sd(i)=1;
    elseif(params_out.RECIST_mets{j}=="PD")
        calc_PD=calc_PD+1;
        vector_pd(i)=1;
    else
        calc_NP=calc_NP+1;
    end

end

calc_PD=calc_PD/npatients_plausible;
calc_CR=calc_CR/npatients_plausible;
calc_PR=calc_PR/npatients_plausible;
calc_SD=calc_SD/npatients_plausible;


nsamples=312;
sample_resp_rates = @(x)mean(datasample(x,nsamples));

CR_CIs = bootci(1000,{sample_resp_rates,vector_cr},'type','per','alpha',.05);
PR_CIs = bootci(1000,{sample_resp_rates,vector_pr},'type','per','alpha',.05);
SD_CIs = bootci(1000,{sample_resp_rates,vector_sd},'type','per','alpha',.05);
PD_CIs = bootci(1000,{sample_resp_rates,vector_pd},'type','per','alpha',.05);


disp('CR')
disp(calc_CR)
disp('CR CIs')
disp(CR_CIs)

disp('PR')
disp(calc_PR)
disp('PR CIs')
disp(PR_CIs)

disp('SD')
disp(calc_SD)
disp('SD CIs')
disp(SD_CIs)

disp('PD')
disp(calc_PD)
disp('PD CIs')
disp(PD_CIs)