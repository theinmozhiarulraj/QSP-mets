% thein: script to calculate the median duration of response (DOR) and median time
% to respond (TTR) in patients with complete/partial response

% reset DOR=0 as NaN to exclude patients with stable or progressive
% disease

DOR_mets = params_out.DOR_mets;
DOR_mets(DOR_mets==0) = NaN;
%DOR_mets(DOR_mets>200) = 800;
median_dor_mets=median(DOR_mets,'omitnan');

disp("Median DOR of mets (months)")
disp(median_dor_mets/30) % convert days to months

% reset TTR=400 as NaN to exclude patients with stable or progressive
% disease

TTR_mets = params_out.TTR_mets;
TTR_mets(TTR_mets==0) = NaN; 
%TTR_mets(TTR_mets==400) = NaN; % not actually needed
median_ttr_mets=median(TTR_mets,'omitnan');

disp("Median TTR of mets (months)")
disp(median_ttr_mets/30) % convert days to months

% confidence intervals for DOR and TTR
nsamples=312;
sample_medians = @(x)median(datasample(x./30,nsamples),'omitnan');
DOR_CIs = bootci(1000,{sample_medians,DOR_mets},'type','per','alpha',.05);
TTR_CIs = bootci(1000,{sample_medians,TTR_mets},'type','per','alpha',.05);

disp("DOR CIs")
disp(DOR_CIs) % convert days to months

disp("TTR CIs")
disp(TTR_CIs) % convert days to months

% number of patients with DOR>=6 months or >=12 months

npatients_DOR_6mo=length(DOR_mets(DOR_mets>=6*30));
npatients_DOR_12mo=length(DOR_mets(DOR_mets>=12*30));

% total number of patients having complete or partial response

total_crpr=length(params_out.RECIST_mets(params_out.RECIST_mets=="CR/PR"));

% fraction of complete and partial responders with DOR>=6 months or >=12
% months

perc_DOR_6mo=npatients_DOR_6mo*100/total_crpr;
perc_DOR_12mo=npatients_DOR_12mo*100/total_crpr;

disp(" Rate of response duration >=6 mo ")
disp(perc_DOR_6mo) % convert days to months

disp(" Rate of response duration >=12 mo ")
disp(perc_DOR_12mo) % convert days to months

writetable(table(DOR_mets),'DOR_allpatients.csv')
writetable(table(TTR_mets),'TTR_allpatients.csv')



