%%
% Calculate subthreshold DELTAS 
%

clear all
close all
clc
%% Load data/ adress folder %% 
[pathname11] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname11]);

load('Subthreshold_juvenline_HVC_deltas_B_G.mat')

contr_mean_trace=mean_call_trace;

contr_list=filelist2;

clearvars -except contr_mean_trace contr_list

load('Subthreshold_juvenline_HVC_deltas_G.mat')

treat_mean_trace=mean_call_trace;

treat_list=filelist2;

clearvars -except contr_mean_trace contr_list treat_mean_trace treat_list

%load('Significant_cell_values_subthreshold_E.mat')



for g=1:length(contr_mean_trace)
    
    temp_avg_contr=mean(contr_mean_trace{g,1}); % demean the trace first
    
    temp_contr_trace=contr_mean_trace{g,1}-temp_avg_contr;
    
    temp_treat_trace=treat_mean_trace{g,1}-temp_avg_contr; % demean it using the control
    
    [max_contr, id_contr]=max(abs(temp_contr_trace));
    
    sign_contr=sign(temp_contr_trace(id_contr));
    
    [max_treat, id_treat]=max(abs(temp_treat_trace));
    
    sign_treat=sign(temp_treat_trace(id_treat));
    
%     figure
%     subplot(2,1,1)
%     plot(temp_contr_trace,'-k')
%     hold on
%     plot([0 length(temp_contr_trace)], [max_contr*sign_contr max_contr*sign_contr], '-b')
%     ylim([min(temp_contr_trace)-2 max(temp_contr_trace)+2])
%     
%     subplot(2,1,2)
%     plot(temp_treat_trace,'-r')
%     hold on
%     plot([0 length(temp_treat_trace)], [max_treat*sign_treat max_treat*sign_treat], '-b')  
%     ylim([min(temp_treat_trace)-2 max(temp_treat_trace)+2])

%     pause
    close all
    
    Delta(g)=(max_treat*sign_treat)-(max_contr*sign_contr);
    
    delta_max_treat(g)=max_treat*sign_treat;
    delta_contr_treat(g)=max_contr*sign_contr;
    
    clear temp_countr_max temp_treat_max idx_contr idx_treat temp_treat_trace temp_contr_trace...
        temp_avg_contr temp_avg_treat max_contr max_treat id_treat id_contr sign_treat sign_contr...
    
end

figure

a=repmat(1,length(Delta),1);

scatter(a,Delta,'mo','Jitter', 'on')
hold on
plot([0 2], [0 0], '--k', 'LineWidth',2)
axis square
box off
xlim([0 2])
ylim([-10 10])

%% do the model analysis here:

%% for E: 
% 
% cell_ID=(1:23).';
% 
% bird_ID=[1,1,1,1,1,1,1,1,2,3,4,5,5,6,6,7,7,8,8,9,9,9,9].';
% 
% bird_ID2=[bird_ID; bird_ID];
% 
% cell_ID2=[cell_ID; cell_ID];
% 
% delta=[delta_contr_treat.'; delta_max_treat.'];
% 
% treatment=[repmat(1,1,length(bird_ID)).'; repmat(2,1,length(bird_ID)).'];
% 
% delta_table=table(bird_ID2, cell_ID2, treatment, delta, 'VariableNames', {'bird_ID', 'cell_ID', 'treatment', 'delta'});
% 
% lme1=fitlme(delta_table,'delta~treatment+(1|bird_ID)+(1|bird_ID:cell_ID)')
% 
% [~,~,stats1]=randomEffects(lme1)
% 
% lme1.Rsquared

%% for G:
% 
% cell_ID=(1:15);
% 
% bird_ID=[1,1,1,1,1,1,1,1,2,3,4,5,4,6,6].';

% bird_ID2=[bird_ID; bird_ID];
% 
% cell_ID2=[cell_ID; cell_ID];
% 
% delta=[delta_contr_treat(cell_idx).'; delta_max_treat(cell_idx).'];
% 
% treatment=[repmat(1,1,length(bird_ID)).'; repmat(2,1,length(bird_ID)).'];
% 
% delta_table=table(bird_ID2, cell_ID2, treatment, delta, 'VariableNames', {'bird_ID', 'cell_ID', 'treatment', 'delta'});
% 
% lme1=fitlme(delta_table,'delta~treatment+(1|bird_ID)+(1|bird_ID:cell_ID)')
% 
% [~,~,stats1]=randomEffects(lme1)
% 
% lme1.Rsquared

%% see whether you can use lme to compare against 0 for E:
% 
% cell_ID_0=(1:23).';
% 
% bird_ID_0=[1,1,1,1,1,1,1,1,2,3,4,5,5,6,6,7,7,8,8,9,9,9,9].';
% 
% 
% zero_0=zeros(23,1);
% 
% delta_table_0=table(bird_ID_0, cell_ID_0, zero_0, Delta.', 'VariableNames', {'bird_ID', 'cell_ID', 'zeros', 'delta'});
% 
% lme_zeros=fitlme(delta_table_0,'delta~0+(1|bird_ID)+(1|bird_ID:cell_ID)')
% 
% [~,~,stats5]=randomEffects(lme_zeros)

%% see whether you can use lme to compare against 0 for G:

cell_ID=(1:15).';

bird_ID=[1,1,1,1,1,1,1,1,2,3,4,5,4,6,6].';

delta_table_0=table(bird_ID, cell_ID, Delta.', 'VariableNames', {'bird_ID', 'cell_ID', 'delta'});

lme_zeros=fitlme(delta_table_0,'delta~0+(1|bird_ID)+(1|bird_ID:cell_ID)')

[~,~,stats5]=randomEffects(lme_zeros)