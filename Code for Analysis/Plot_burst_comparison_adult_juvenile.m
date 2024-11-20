clear all 
close all
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);

%% load data:

load('Adult_HVC_data_during_singing_burst_isi_analysis8.mat')

adu_filelist=filelist2;

a_burst_spikes=cell2mat(cell_spikes_per_burst);

a_spikes_per_burst=cell_spikes_per_burst;

a_cell_burst_time=cell_burst_time;

a_max_bursts=nr_bursts_max;

a_burst_duration=cell_burst_duration;

a_burst_inst_fr=cell_burst_inst_firing_rate;

a_burst_perc_on=cell_burst_percent_onset;

a_burst_per_dur=cell_burst_percent_duration;

a_single_spikes=cell_single_spikes;

a_burst_isi=cell_burst_isi_duration;

a_nr_spikes_trial_mean=nr_spikes_trial_mean;


clearvars -except a_burst_spikes a_max_bursts a_cell_burst_time a_spikes_per_burst a_burst_duration a_burst_inst_fr a_burst_perc_on...
    a_burst_per_dur a_single_spikes a_cell_amp a_freq_cells a_burst_isi adu_filelist a_nr_spikes_trial_mean

load('Adult_HVC_data_psp_duration_analysis2.mat')


a_cell_amp=average_amplitude_per_cell;

a_freq_cells=frequency_per_cell;

a_psp_duration=psp_durations_cells;

a_psp_duration_corrected=psp_durations_cells_corrected;
a_psp_amplitude_cells_corrected=psp_amplitude_cells_corrected;

clearvars -except a_burst_spikes a_max_bursts a_cell_burst_time a_spikes_per_burst a_burst_duration a_burst_inst_fr a_burst_perc_on...
    a_burst_per_dur a_single_spikes a_cell_amp a_freq_cells a_burst_isi a_psp_duration adu_filelist a_psp_duration_corrected...
    a_psp_amplitude_cells_corrected a_nr_spikes_trial_mean



load('Juvenile_HVC_data_psp_duration_analysis3.mat')

j_cell_amp=average_amplitude_per_cell;

j_freq_cells=frequency_per_cell;

j_psp_duration=psp_durations_cells;

j_psp_duration_corrected=psp_durations_cells_corrected;
j_psp_amplitude_cells_corrected=psp_amplitude_cells_corrected;

clearvars -except a_burst_spikes a_max_bursts j_burst_spikes j_max_bursts a_cell_burst_time j_cell_burst_time a_spikes_per_burst j_spikes_per_burst...
    j_burst_variance2  a_burst_duration a_burst_inst_fr j_burst_inst_fr j_burst_per_dur j_burst_perc_on a_burst_per_dur a_burst_perc_on...
    a_single_spikes j_single_spikes a_cell_amp a_freq_cells j_cell_amp j_freq_cells a_burst_isi a_psp_duration j_psp_duration adu_filelist...
    a_psp_duration_corrected j_psp_duration_corrected a_psp_amplitude_cells_corrected j_psp_amplitude_cells_corrected a_nr_spikes_trial_mean


load('Juvenile_HVC_data_during_singing_burst_isi_analysis11.mat')

j_burst_spikes=cell2mat(cell_spikes_per_burst);

j_spikes_per_burst=cell_spikes_per_burst;

j_cell_burst_time=cell_burst_time;

j_max_bursts=nr_bursts_max;

%j_burst_variance2=burst_variance;

juv_filelist=filelist2;

j_burst_inst_fr=cell_burst_inst_firing_rate;

j_burst_perc_on=cell_burst_percent_onset;

j_burst_per_dur=cell_burst_percent_duration;

j_single_spikes=cell_single_spikes;

j_burst_duration=cell_burst_duration;

j_burst_isi=cell_burst_isi_duration;

j_nr_spikes_trial_mean=nr_spikes_trial_mean;

clearvars -except a_burst_spikes a_max_bursts j_burst_spikes j_max_bursts a_cell_burst_time j_cell_burst_time a_spikes_per_burst j_spikes_per_burst...
    j_burst_variance2  a_burst_duration a_burst_inst_fr j_burst_inst_fr j_burst_per_dur j_burst_perc_on a_burst_per_dur a_burst_perc_on...
    a_single_spikes j_single_spikes a_cell_amp a_freq_cells j_cell_amp j_freq_cells j_burst_duration a_burst_isi j_burst_isi a_psp_duration j_psp_duration adu_filelist juv_filelist...
    a_psp_duration_corrected j_psp_duration_corrected a_psp_amplitude_cells_corrected j_psp_amplitude_cells_corrected a_nr_spikes_trial_mean...
    j_nr_spikes_trial_mean

load('MD_birds_syllable_ISI.mat')

ISI_filelist=filelist;

clear filelist

%%

% first make a list of bird names per cell and cell IDs:

for z = 1:length(adu_filelist)
    
    
name = char(adu_filelist(z,1).name);



if z>32
    
    adu_bird{1,z}=name(1:6);
    
    
else
    
    adu_bird{1,z}=name(1:4);
    
end

clear name 
end

% then same for juveniles:

for z = 1:length(juv_filelist)
    
    
name = char(juv_filelist(z,1).name);



if z>7
    
    juv_bird{1,z}=name(10:13);
    
    
else
    
    juv_bird{1,z}=name(12:14);
    
end

clear name 
end

a_cell_ID=[1:length(adu_bird)];

j_cell_ID=[1:length(juv_bird)];

str1='juvenile';
str2='adult';

%%
% for g=1:length(j_burst_spike_variance)
% 
% j_burst_var(g)=std(j_burst_spike_variance{1,g})
% 
% end
% 
% for g=1:length(a_burst_spike_variance)
% 
% a_burst_var(g)=std(a_burst_spike_variance{1,g})
% 
% end
% 
% max(a_burst_spikes)
% max(j_burst_spikes)
% 
% edges1=0:1:max(a_burst_spikes);
% 
% a_bin_idx1=discretize(a_burst_spikes,edges1,'IncludedEdge','right');
% j_bin_idx1=discretize(j_burst_spikes,edges1,'IncludedEdge','right');
% 
% a_groups1=unique(a_bin_idx1);
% 
% for t=1:length(a_groups1)
%    
%     a1_howmany=find(a_groups1(t)==a_bin_idx1);
%     a1_percent(t)=(length(a1_howmany)/length(a_burst_spikes))*100;
%     
%     clear howmany
% end
% 
% j_groups1=unique(j_bin_idx1);
% 
% for t=1:length(j_groups1)
%    
%     j1_howmany=find(j_groups1(t)==j_bin_idx1);
%     j1_percent(t)=(length(j1_howmany)/length(j_burst_spikes))*100;
%     
%     clear howmany
% end
% 
% figure
% 
% bar(a_groups1,a1_percent,'k')
% hold on
% bar(j_groups1,j1_percent,'r')
% ylabel('% Cells')
% xlabel('spikes per burst')
% axis square
% box off


%% plot max bursts per motif

max(a_max_bursts)
max(j_max_bursts)

edges=-0.5:1:max(a_max_bursts+0.5);

a_bin_idx=discretize(a_max_bursts,edges,'IncludedEdge','right');
j_bin_idx=discretize(j_max_bursts,edges,'IncludedEdge','right');

a_groups=unique(a_bin_idx);


for t=1:length(a_groups)
   
    a_howmany=find(a_groups(t)==a_bin_idx);
    a_percent(t)=(length(a_howmany)/length(a_max_bursts))*100;
    
    clear howmany
end

j_groups=unique(j_bin_idx);

for t=1:length(j_groups)
   
    j_howmany=find(j_groups(t)==j_bin_idx);
    j_percent(t)=(length(j_howmany)/length(j_max_bursts))*100;
    
    clear howmany
end

% make sure to plot each edge:

j_mat=zeros(1,length(edges));

j_mat(j_groups)=j_percent;

j_mat2=[0, j_mat, 0];

a_mat=zeros(1,length(edges));

a_mat(a_groups)=a_percent;

a_mat2=[0, a_mat, 0];

% for stairs add zeros at the beginning and end:

edges2=[-0.5,edges,max(edges)+0.5];

figure

stairs(edges2,a_mat2,'k')
hold on
stairs(edges2,j_mat2,'r')
hold on
plot([mean(j_max_bursts) mean(j_max_bursts)], [50 50], 'rv')
hold on
plot([mean(a_max_bursts) mean(a_max_bursts)], [50 50], 'kv')
ylabel('% Cells')
xlabel('Max Bursts')
xlim([-0.5 5])
axis square
box off

% do a kolomogrov smirnov test for comparing distributions:

[h,p]=kstest2(a_mat, j_mat); % for continuous data
[tbl, chi2, p]=crosstab(a_mat,j_mat); % for discrete data (like histograms)


clear tbl chi2 p

clear edges2 a_mat2 j_mat2 j_mat a_mat edges h p


%% plot spikes per burst


max(a_burst_spikes)
max(j_burst_spikes)

% clean from NaNs:

ida=find(isnan(a_burst_spikes)==1);
a_burst_spikes(ida)=[];

idu=find(isnan(j_burst_spikes)==1);
j_burst_spikes(idu)=[];

edges_sp=0.5:1:max(a_burst_spikes)+0.5;

a_bin_idx_sp=discretize(a_burst_spikes,edges_sp,'IncludedEdge','right');
j_bin_idx_sp=discretize(j_burst_spikes,edges_sp,'IncludedEdge','right');

a_groups_sp=unique(a_bin_idx_sp);


for t=1:length(a_groups_sp)
   
    a_howmany_sp=find(a_groups_sp(t)==a_bin_idx_sp);
    a_percent_sp(t)=(length(a_howmany_sp)/length(a_burst_spikes))*100;
    
    clear howmany
end

j_groups_sp=unique(j_bin_idx_sp);

for t=1:length(j_groups_sp)
   
    j_howmany_sp=find(j_groups_sp(t)==j_bin_idx_sp);
    j_percent_sp(t)=(length(j_howmany_sp)/length(j_burst_spikes))*100;
    
    clear howmany
end

% make sure to plot each edge:

j_mate=zeros(1,length(edges_sp));

j_mate(j_groups_sp)=j_percent_sp;

j_mate2=[0, j_mate, 0];

a_mate=zeros(1,length(edges_sp));

a_mate(a_groups_sp)=a_percent_sp;

a_mate2=[0, a_mate, 0];

% for stairs add zeros at the beginning and end:

edges2=[0,edges_sp,max(edges_sp)+0.5];

figure

stairs(edges2,a_mate2,'k')
hold on
stairs(edges2,j_mate2,'r')
hold on
plot([mean(j_burst_spikes) mean(j_burst_spikes)], [40 40], 'rv')
hold on
plot([mean(a_burst_spikes) mean(a_burst_spikes)], [40 40], 'kv')
ylabel('% Bursts')
xlabel('Spikes per burst')
xlim([1 12])
axis square
box off

clear edges2 a_mate2 j_mate2 j_mate a_mate edges




%% find corresponding bursts:
 %burst is considered the same within 20 ms boundary

% adults:
 
for g=1:length(a_cell_burst_time)
    
    idx=0;
   
    who=(a_cell_burst_time{1,g});
    
    indeces=NaN(1,length(who));
    
    for h=1:length(who)
        
            idx=idx+1;    
       
        for f=(h+1):length(who)
        
            bet=abs(who(f)-who(h));
        
                if bet<20
            
                    if isnan(indeces(f))==1
                    
                    indeces(f)=idx;
                    
                    end
                    
                    if isnan(indeces(h))==1
                    
                    indeces(h)=idx;
                    
                    end
                    
                end
        
        end
        
    end
    
    indeces2=indeces;
    
    frog=find(isnan(indeces2)==1); % find the NaN values (unique bursts) and mark them with an idx:
    
    for a=1:length(frog)
       
       idx=idx+1;
       indeces2(frog(a))=idx;
        
    end
    
    a_burst_time_idx{1,g}=indeces;
    a_burst_time_idx_all{1,g}=indeces2;
    
    clear indeces who indeces2 frog
    
end

% juveniles:

for g=1:length(j_cell_burst_time)
    
    idx=0;
   
    who=(j_cell_burst_time{1,g});
    
    indeces=NaN(1,length(who));
    
    for h=1:length(who)
        
            idx=idx+1;    
       
        for f=(h+1):length(who)
        
            bet=abs(who(f)-who(h));
        
                if bet<20
            
                    if isnan(indeces(f))==1
                    
                    indeces(f)=idx;
                    
                    end
                    
                    if isnan(indeces(h))==1
                    
                    indeces(h)=idx;
                    
                    end
                    
                end
        
        end
        
    end
    
        
    indeces2=indeces;
    
    frog=find(isnan(indeces2)==1); % find the NaN values (unique bursts) and mark them with an idx:
    
    for a=1:length(frog)
       
       idx=idx+1;
       indeces2(frog(a))=idx;
        
    end
    
    j_burst_time_idx{1,g}=indeces;
    j_burst_time_idx_all{1,g}=indeces2;
    
    clear indeces who indeces2 frog
    
end

%% calculate burst duration:

% adults:

c1=0;

for t=1:length(a_burst_duration)
    
        temp_spikes=a_spikes_per_burst{1,t};
        temp_dur=a_burst_duration{1,t};
        temp_idx=a_burst_time_idx_all{1,t};
        str=unique(a_burst_time_idx_all{1,t});
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];

        if isempty(str)==0
        
                for z=1:length(str)
                    
                    c1=c1+1;
            
                    a_all_burst_bird{1,c1}=adu_bird{1,t};
                    a_all_burst_cell_ID(c1)=a_cell_ID(t);
                    

                    id=find(str(z)==temp_idx);
                    
                    mean_spikes_per_burst(z)=mean(temp_spikes(id))
                    
                    burst_duration_mean(z)=mean(temp_dur(id));
                    

                    clear id
                end

            a_burst_mean_duration{1,t}=burst_duration_mean;
            a_mean_spikes_per_burst{1,t}=mean_spikes_per_burst;
    
        else
            
            a_burst_mean_duration{1,t}=NaN;
            a_mean_spikes_per_burst{1,t}=NaN;
            
        end
    clear temp_dur temp_idx str u burst_duration_mean temp_spikes burst_duration_mean mean_spikes_per_burst
end


a_burst_mean_mat=cell2mat(a_burst_mean_duration);
a_mean_spikes_per_burst_mat=cell2mat(a_mean_spikes_per_burst);

no=find(isnan(a_burst_mean_mat)==1);

a_burst_mean_mat(no)=[];
a_mean_spikes_per_burst_mat(no)=[];

clear no

% juveniles: 

c1=0;

for t=1:length(j_burst_duration)
        
        temp_spikes=j_spikes_per_burst{1,t};
        temp_dur=j_burst_duration{1,t};
        temp_idx=j_burst_time_idx_all{1,t};
        str=unique(j_burst_time_idx_all{1,t});
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];

        if isempty(str)==0
        
                for z=1:length(str)
                    
                    c1=c1+1;
            
                    j_all_burst_bird{1,c1}=juv_bird{1,t};
                    j_all_burst_cell_ID(c1)=j_cell_ID(t);

                    id=find(str(z)==temp_idx);
                    
                    mean_spikes_per_burst(z)=mean(temp_spikes(id))
                    
                    burst_duration_mean(z)=mean(temp_dur(id));
                    

                    clear id
                end

            j_burst_mean_duration{1,t}=burst_duration_mean;
            j_mean_spikes_per_burst{1,t}=mean_spikes_per_burst;
    
        else
            
            j_burst_mean_duration{1,t}=NaN;
            j_mean_spikes_per_burst{1,t}=NaN;
            
        end
    clear temp_dur temp_idx str u burst_duration_mean temp_dur temp_idx str u...
        burst_duration_mean temp_spikes burst_duration_mean mean_spikes_per_burst
end


j_burst_mean_mat=cell2mat(j_burst_mean_duration);
j_mean_spikes_per_burst_mat=cell2mat(j_mean_spikes_per_burst);

no=find(isnan(j_burst_mean_mat)==1);

j_burst_mean_mat(no)=[];
j_mean_spikes_per_burst_mat(no)=[];

clear no

% test1=1/j_burst_mean_mat(:);
% test2=1/a_burst_mean_mat(:);


%% calculate variance per burst
% 
% a_burst_variance=[];
% 
% for t=1:length(a_spikes_per_burst)
%     
%     counter=0;
% 
%         temp_spikes=a_spikes_per_burst{1,t};
%         temp_idx=a_burst_time_idx{1,t};
%         str=unique(a_burst_time_idx{1,t});
%         sanity_check{1,t}=str;
% 
%         for z=1:length(str)
%    
%             id=find(str(z)==temp_idx);
%         
%             if length(id)>1
%             
%                 counter=counter+1;
%                 
%                 temp_var=var(temp_spikes(id));
%                 a_burst_variance{1,t}(counter)=temp_var;
%             
%             end
%         
%             clear id temp_var
%         end
%     
%     clear temp_spikes temp_idx str
% end
% 
% a_burst_var_mat=cell2mat(a_burst_variance);
% 
% j_burst_variance=[];
% 
% for t=1:length(j_spikes_per_burst)
%     
%     counter=0;
% 
%         temp_spikes=j_spikes_per_burst{1,t};
%         temp_idx=j_burst_time_idx{1,t};
%         str=unique(j_burst_time_idx{1,t});
%         sanity_check{1,t}=str;
% 
%         for z=1:length(str)
%    
%             id=find(str(z)==temp_idx);
%         
%             if length(id)>1
%             
%                 counter=counter+1;
%                 
%                 temp_var=var(temp_spikes(id));
%                 j_burst_variance{1,t}(counter)=temp_var;
%             
%             end
%         
%             clear id temp_var
%         end
%     
%     clear temp_spikes temp_idx str
% end
% 
% j_burst_var_mat=cell2mat(j_burst_variance);

%% plot burst variance:
% 
% max(a_burst_var_mat)
% max(j_burst_variance2)
% 
% edges2=0:0.1:max(j_burst_variance2)+1;
% 
% a_bin_idx2=discretize(a_burst_var_mat,edges2,'IncludedEdge','right');
% j_bin_idx2=discretize(j_burst_var_mat,edges2,'IncludedEdge','right');
% 
% a_groups2=unique(a_bin_idx2);
% 
% 
% for t=1:length(a_groups2)
%    
%     a_howmany2=find(a_groups2(t)==a_bin_idx2);
%     a_percent2(t)=(length(a_howmany2)/length(a_burst_var_mat))*100;
%     
%     clear howmany
% end
% 
% j_groups2=unique(j_bin_idx2);
% 
% for t=1:length(j_groups2)
%    
%     j_howmany2=find(j_groups2(t)==j_bin_idx2);
%     j_percent2(t)=(length(j_howmany2)/length(j_burst_var_mat))*100;
%     
%     clear howmany
% end
% 
% 
% figure
% 
% bar(edges2(a_groups2),a_percent2,'k')
% hold on
% bar(edges2(j_groups2),j_percent2,'r')
% ylabel('% bursts')
% xlabel('burst variance')
% axis square
% box off

% do 2 subplots with hist:

% figure
% subplot(1,2,1)
% hist(a_burst_var_mat)
% ylim([0 length(a_burst_var_mat)])
% ytix = get(gca, 'YTick')
% set(gca, 'YTick',ytix, 'YTickLabel',ytix/length(a_burst_var_mat)*100) % change y axis to percentage here
% 
% 
% subplot(1,2,2)
% hist(j_burst_var_mat)
% ylim([0 length(j_burst_var_mat)])
% ytix2 = get(gca, 'YTick')
% set(gca, 'YTick',ytix2, 'YTickLabel',ytix2/length(j_burst_var_mat)*100) % change y axis to percentage here

%% calculate burst jitter:
% modus (mode) - most frequent value in an array

% adults:

% a_burst_mode=[];

c1=0;

for t=1:length(a_spikes_per_burst)


        temp_spikes=a_spikes_per_burst{1,t};
        temp_idx=a_burst_time_idx{1,t};
        str=unique(a_burst_time_idx{1,t});
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];

        if isempty(str)==0
        

            
            
                for z=1:length(str)

                    id=find(str(z)==temp_idx);

                    c1=c1+1;
            
                    a_repeated_burst_bird{1,c1}=adu_bird{1,t};
                    a_repeated_burst_cell_ID(c1)=a_cell_ID(t);
                    
                    
                    common_val=mode(temp_spikes(id));
                    jitter=abs(temp_spikes(id)-common_val);
                    a_burst_jitter{1,t}{1,z}=jitter;
                    temp_mean_jitter(z)=mean(jitter);
                    temp_max_jitter(z)=max(jitter);
                    
%                     a_burst_mode=[a_burst_mode, common_val];

                    clear id temp_var common_val jitter
                end

            a_mean_burst_jitters{1,t}=temp_mean_jitter;
            a_max_burst_jitters{1,t}=temp_max_jitter;

           
    
        else
            
            a_mean_burst_jitters{1,t}=NaN;
            a_max_burst_jitters{1,t}=NaN;

            
        end
    clear temp_spikes temp_idx str temp_mean_jitter temp_max_jitter mean_spikes_per_burst
end


a_jitters_mat=cell2mat(a_mean_burst_jitters);

idm=find(isnan(a_jitters_mat)==1);

a_jitters_mat(idm)=[]; % get rid of the NaN values

clear idm idn ida idu

a_max_jitters_mat=cell2mat(a_max_burst_jitters);

idn=find(isnan(a_max_jitters_mat)==1);

a_max_jitters_mat(idn)=[]; % get rid of the NaN values

clear idn


% juveniles:

% j_burst_mode=[];

c1=0;

for t=1:length(j_spikes_per_burst)


        temp_spikes=j_spikes_per_burst{1,t};
        temp_idx=j_burst_time_idx{1,t};
        str=unique(j_burst_time_idx{1,t});
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];

        if isempty(str)==0
            

        
                for z=1:length(str)

                    id=find(str(z)==temp_idx);
                    
                    c1=c1+1;
            
                    j_repeated_burst_bird{1,c1}=juv_bird{1,t};
                    j_repeated_burst_cell_ID(c1)=j_cell_ID(t);

                    common_val=mode(temp_spikes(id));
                    jitter=abs(temp_spikes(id)-common_val);
                    j_burst_jitter{1,t}{1,z}=jitter;
                    temp_mean_jitter(z)=mean(jitter);
                    temp_max_jitter(z)=max(jitter);
%                     j_burst_mode=[j_burst_mode, common_val];


                    clear id temp_var common_val jitter
                end

            j_mean_burst_jitters{1,t}=temp_mean_jitter;
            j_max_burst_jitters{1,t}=temp_max_jitter;
    
        else
            
            j_mean_burst_jitters{1,t}=NaN;
            j_max_burst_jitters{1,t}=NaN;
            
        end
    clear temp_spikes temp_idx str temp_mean_jitter temp_max_jitter
end


j_jitters_mat=cell2mat(j_mean_burst_jitters);

idn=find(isnan(j_jitters_mat)==1);

j_jitters_mat(idn)=[]; % get rid of the NaN values

clear idm idn ida idu

j_max_jitters_mat=cell2mat(j_max_burst_jitters);

idn=find(isnan(j_max_jitters_mat)==1);

j_max_jitters_mat(idn)=[]; % get rid of the NaN values

clear idn

figure

subplot(1,2,1)

hist(a_jitters_mat)

subplot(1,2,2)

hist(j_jitters_mat)




figure

subplot(1,2,1)

hist(a_max_jitters_mat)

subplot(1,2,2)

hist(j_max_jitters_mat)

%% plot burst jitter:

clear edges edges1 deges2 edges_sp a_groups j_groups a_howmany a_percent j_groups j_percent j_howmany...
    a_groups_stairs a_percent_stairs j_groups_stairs j_percent_stairs

max(a_jitters_mat)
max(j_jitters_mat)

edges=0:0.1:max(j_jitters_mat)+0.1;

a_bin=discretize(a_jitters_mat,edges,'IncludedEdge','right');
j_bin=discretize(j_jitters_mat,edges,'IncludedEdge','right');

a_groups=unique(a_bin);


for t=1:length(a_groups)
   
    a_howmany=find(a_groups(t)==a_bin);
    a_percent(t)=(length(a_howmany)/length(a_jitters_mat))*100;
    
    clear howmany
end

j_groups=unique(j_bin);

for t=1:length(j_groups)
   
    j_howmany=find(j_groups(t)==j_bin);
    j_percent(t)=(length(j_howmany)/length(j_jitters_mat))*100;
    
    clear howmany
end


% make sure to plot each edge:

j_mate=zeros(1,length(edges));

j_mate(j_groups)=j_percent;

j_mate2=[0, j_mate, 0];

a_mate=zeros(1,length(edges));

a_mate(a_groups)=a_percent;

a_mate2=[0, a_mate, 0];

% for stairs add zeros at the beginning and end:

edges2=[0,edges,max(edges)+0.1];

figure

stairs(edges2,a_mate2,'k')
hold on
stairs(edges2,j_mate2,'r')
hold on
plot([mean(j_jitters_mat) mean(j_jitters_mat)], [60 60], 'rv')
hold on
plot([mean(a_jitters_mat) mean(a_jitters_mat)], [60 60], 'kv')
ylabel('% Bursts')
xlabel('Mean Jitter (Spikes per burst)')
xlim([-0.2 max(edges)])
axis square
box off


[tbl, chi2, p]=crosstab(a_mate, j_mate)

clear tbl chi2 p

%% calculate instantaneous firing rate:

% adults:

for t=1:length(a_burst_inst_fr)

        temp_inst_fr=a_burst_inst_fr{1,t};
        temp_idx=a_burst_time_idx_all{1,t};
        str=unique(a_burst_time_idx_all{1,t});
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];
        burst_isi_inv=[];

        if isempty(str)==0
        
                for z=1:length(str)

                    id=find(str(z)==temp_idx);
                    
                    burst_inst_fr_mean(z)=mean(temp_inst_fr(id));
                    

                    clear id
                end

            a_inst_fr_burst{1,t}=burst_inst_fr_mean;
    
        else
            
            a_inst_fr_burst{1,t}=NaN;
            
        end
    clear temp_nr_spikes temp_inst_fr burst_isi_inv mean_burst_fr temp_idx str burst_isi_inv_mean burst_inst_fr_mean
end


a_mean_burst_fr_mat=cell2mat(a_inst_fr_burst);


no=find(isnan(a_mean_burst_fr_mat)==1);

a_mean_burst_fr_mat(no)=[];

% juveniles:

for t=1:length(j_burst_inst_fr)

        temp_inst_fr=j_burst_inst_fr{1,t};
        temp_idx=j_burst_time_idx_all{1,t};
        str=unique(j_burst_time_idx_all{1,t});
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];
        burst_isi_inv=[];

        if isempty(str)==0
        
                for z=1:length(str)

                    id=find(str(z)==temp_idx);
                    
                    burst_inst_fr_mean(z)=mean(temp_inst_fr(id));

                    clear id
                end

            j_inst_fr_burst{1,t}=burst_inst_fr_mean;
    
        else
            
            j_inst_fr_burst{1,t}=NaN;
            
        end
    clear temp_nr_spikes temp_inst_fr burst_isi_inv mean_burst_fr temp_idx str burst_isi_inv_mean burst_inst_fr_mean
end


j_mean_burst_fr_mat=cell2mat(j_inst_fr_burst);


no=find(isnan(j_mean_burst_fr_mat)==1);

j_mean_burst_fr_mat(no)=[];

juv_scatter=repmat(1,length(j_mean_burst_fr_mat),1);

adt_scatter=repmat(2,length(a_mean_burst_fr_mat),1);

figure

subplot(1,2,1)
scatter(juv_scatter, j_mean_burst_fr_mat,'ro','Jitter','on')
hold on
plot([0.8 1.2],[mean(j_mean_burst_fr_mat) mean(j_mean_burst_fr_mat)], '-k')
hold on
scatter(adt_scatter, a_mean_burst_fr_mat,'bo','Jitter','on')
hold on
plot([1.8 2.2],[mean(a_mean_burst_fr_mat) mean(a_mean_burst_fr_mat)], '-k')
hold on
ylim([100 800])
xlim([0 3])
ylabel('Instantaneous Firing Rate (Hz)')
axis square
box off


figure

subplot(1,2,1)
hist(a_mean_burst_fr_mat)
xlim([0 800])

subplot(1,2,2)
hist(j_mean_burst_fr_mat)
xlim([0 800])


max(a_mean_burst_fr_mat)
max(j_mean_burst_fr_mat)

clear a_bin_idx j_bin_idx a_groups a_howmany a_percent j_groups j_howmany j_percent edges

edges=0:10:max(a_mean_burst_fr_mat+10);

a_bin_idx=discretize(a_mean_burst_fr_mat,edges,'IncludedEdge','right');
j_bin_idx=discretize(j_mean_burst_fr_mat,edges,'IncludedEdge','right');

a_groups=unique(a_bin_idx);

a_groups(isnan(a_groups))=[]; % remove the NaNs


for t=1:length(a_groups)
   
    a_howmany=find(a_groups(t)==a_bin_idx);
    a_percent(t)=(length(a_howmany)/length(a_mean_burst_fr_mat))*100;
    
    clear howmany
end

j_groups=unique(j_bin_idx);

j_groups(isnan(j_groups))=[]; % remove the NaNs

for t=1:length(j_groups)
   
    j_howmany=find(j_groups(t)==j_bin_idx);
    j_percent(t)=(length(j_howmany)/length(j_mean_burst_fr_mat))*100;
    
    clear howmany
end

% make sure to plot each edge:

j_mat=zeros(1,length(edges));

j_mat(j_groups)=j_percent;

j_mat2=[0, j_mat, 0];

a_mat=zeros(1,length(edges));

a_mat(a_groups)=a_percent;

a_mat2=[0, a_mat, 0];

% for stairs add zeros at the beginning and end:

edges2=[-0.5,edges,max(edges)+0.5];

figure

stairs(edges2,a_mat2,'k')
hold on
stairs(edges2,j_mat2,'r')
hold on
plot([mean(j_mean_burst_fr_mat) mean(j_mean_burst_fr_mat)], [50 50], 'rv')
hold on
plot([mean(a_mean_burst_fr_mat) mean(a_mean_burst_fr_mat)], [50 50], 'kv')
ylabel('% Bursts')
xlabel('Instantaneous firing rate')
xlim([-0.5 800])
ylim([0 20])
axis square
box off

% do a kolomogrov smirnov test for comparing distributions:

%[h,p]=kstest2(a_mat, j_mat); % for continuous data
[tbl, chi2, p]=crosstab(a_mat,j_mat) % for discrete data (like histograms, all binned data)


clear tbl chi2 p

clear edges2 a_mat2 j_mat2 j_mat a_mat edges h p

clear a_bin_idx j_bin_idx a_groups a_howmany a_percent j_groups j_howmany j_percent

figure

scatter(repmat(1,length(j_mean_burst_fr_mat),1), j_mean_burst_fr_mat, 'ro', 'Jitter', 'on')
hold on
scatter(repmat(2,length(a_mean_burst_fr_mat),1), a_mean_burst_fr_mat, 'ko', 'Jitter', 'on')
hold on
plot([0.8 1.2], [mean(j_mean_burst_fr_mat) mean(j_mean_burst_fr_mat)], '-k')
plot([1.8 2.2], [mean(a_mean_burst_fr_mat) mean(a_mean_burst_fr_mat)], '-k')
xlim([0 3])
axis square
box off
ylabel('Burst instantaneous firing rate')
%% sort and plot bursts in percentage:


for t=1:length(j_burst_per_dur)
    
        temp_per_dur=j_burst_per_dur{1,t};
        temp_per_on=j_burst_perc_on{1,t};

        temp_idx=j_burst_time_idx_all{1,t};
        str=unique(temp_idx);
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];

        if isempty(str)==0
        
                for z=1:length(str)

                    id=find(str(z)==temp_idx);
                    
                    burst_perc_on_mean(z)=mean(temp_per_on(id));
                    burst_perc_dur_mean(z)=mean(temp_per_dur(id));
                    

                    clear id temp_mode temp_jitter
                end

            j_burst_mean_perc_on{1,t}=burst_perc_on_mean;
            j_burst_mean_perc_dur{1,t}=burst_perc_dur_mean;
    
        else
            continue
%             j_burst_mean_perc_on{1,t}=NaN;
%             j_burst_mean_perc_dur{1,t}=NaN;
            
        end
    clear temp_per_dur temp_per_on burst_perc_on_mean burst_perc_dur_mean temp_idx str burst_onset_jitter temp_ms_on
end

sorting=[];

for t=1:length(j_burst_mean_perc_on)

    if isempty(j_burst_mean_perc_on{1,t})==0; % check that the cell is not empty
        
        temp_sort=sort(j_burst_mean_perc_on{1,t});
    
        sorting=[sorting, temp_sort(1)];
        
        clear temp_sort
        
    else
        
        sorting=[sorting,NaN];

    end

end

[sorted,j_idx]=sort(sorting,'descend');

figure

subplot(1,2,1)

for t=1:length(j_idx)
       
    temp_rep=repmat(t,1,length(j_burst_mean_perc_on{1,j_idx(t)}));
       
       plot(j_burst_mean_perc_on{1,j_idx(t)},temp_rep,'mo')
       hold on
       
      clear temp_rep
      
end

% adults:

for t=1:length(a_burst_per_dur)
    
        temp_per_dur=a_burst_per_dur{1,t};
        temp_per_on=a_burst_perc_on{1,t};

        temp_idx=a_burst_time_idx_all{1,t};
        str=unique(temp_idx);
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];

        if isempty(str)==0
        
                for z=1:length(str)

                    id=find(str(z)==temp_idx);
                    
                    burst_perc_on_mean(z)=mean(temp_per_on(id));
                    burst_perc_dur_mean(z)=mean(temp_per_dur(id));
                    


                    clear id temp_mode temp_jitter

                end
                

            a_burst_mean_perc_on{1,t}=burst_perc_on_mean;
            a_burst_mean_perc_dur{1,t}=burst_perc_dur_mean;
    
        else
            continue
%             j_burst_mean_perc_on{1,t}=NaN;
%             j_burst_mean_perc_dur{1,t}=NaN;
            
        end
    clear temp_per_dur temp_per_on burst_perc_on_mean burst_perc_dur_mean temp_idx str burst_onset_jitter temp_ms_on
end

a_sorting=[];

for t=1:length(a_burst_mean_perc_on)

    if isempty(a_burst_mean_perc_on{1,t})==0; % check that the cell is not empty
        
        temp_sort=sort(a_burst_mean_perc_on{1,t});
    
        a_sorting=[a_sorting, temp_sort(1)];
        
        clear temp_sort
        
    else
        
        a_sorting=[a_sorting,NaN];

    end

end

[a_sorted,a_idx]=sort(a_sorting,'descend');


subplot(1,2,2)

for t=1:length(a_idx)
       
    temp_rep=repmat(t,1,length(a_burst_mean_perc_on{1,a_idx(t)}));
       
       plot(a_burst_mean_perc_on{1,a_idx(t)},temp_rep,'mo')
       hold on
       
      clear temp_rep
      
end

%% plot burst onsets and offsets:

% juveniles:

figure



j_mega_matrix=[];
j_count=0;


for h=1:length(j_idx)
    
    if isempty(j_burst_mean_perc_on{1,j_idx(h)})==0; % check that the cell is not empty
    
        temp_burst_on=round(j_burst_mean_perc_on{1,j_idx(h)});
        
        temp_small_idx=find(j_burst_mean_perc_dur{1,j_idx(h)}<1);
        
        temp_burst_dur=round(j_burst_mean_perc_dur{1,j_idx(h)});
        
            if isempty(temp_small_idx)==0 % make sure that every burst has at least 1% duration (lowest unit in plot is 1%)

                temp_burst_dur(temp_small_idx)=1;

            end
        
        temp_burst_off=round(temp_burst_on+temp_burst_dur);
        
        j_count=j_count+1;
        
        j_vector=zeros(1, 100); % one vector per cell

            for t=1:length(temp_burst_on)

                temp_fill=ones(1,temp_burst_dur(t));
                j_vector(temp_burst_on(t):temp_burst_off(t)-1)=temp_fill;
                
                temp_rep=repmat(j_count,1,sum(temp_fill));
                
                
                subplot(2,2,1)
                plot([temp_burst_on(t) temp_burst_off(t)-1], [j_count j_count], '-or', 'LineWidth', 5)
                xlim([0 100])
                ylim([0 9])
                hold on

                clear temp_fill temp_rep
            end

        j_mega_matrix=[j_mega_matrix; j_vector];
        
        clear temp_burst_on temp_burst_dur temp_burst_off j_vector
        
    end
    
    clear temp_burst_on temp_burst_dur temp_burst_off j_vector
end

xlabel('Motif duration (%)')
ylabel('Cell ID')

j_summed=sum(j_mega_matrix);
j_summed_prob=smooth(j_summed./j_count); % smoothed using 5 point average


subplot(2,2,3)

plot(j_summed_prob, '-r', 'LineWidth', 2)
ylabel('Burst probability')
xlabel('Motif duration (%)')
xlim([0 100])
ylim([0 0.13])
% check integral function for summation to get 1

% adults:

a_mega_matrix=[];
a_count=0;


for h=1:length(a_idx)
    
    if isempty(a_burst_mean_perc_on{1,a_idx(h)})==0; % check that the cell is not empty
    
        temp_burst_on=round(a_burst_mean_perc_on{1,a_idx(h)});
        temp_burst_dur=round(a_burst_mean_perc_dur{1,a_idx(h)});
        
        temp_small_idx=find(a_burst_mean_perc_dur{1,a_idx(h)}<1);
        
            if isempty(temp_small_idx)==0 % make sure that every burst has at least 1% duration (lowest unit in plot is 1%)

                temp_burst_dur(temp_small_idx)=1;

            end
        
        temp_burst_off=round(temp_burst_on+temp_burst_dur);
        
        a_count=a_count+1;
        
        a_vector=zeros(1, 100); % one vector per cell

            for t=1:length(temp_burst_on)

                temp_fill=ones(1,temp_burst_dur(t));
                a_vector(temp_burst_on(t):temp_burst_off(t)-1)=temp_fill;
                
  
                subplot(2,2,2)
                plot([temp_burst_on(t) temp_burst_off(t)-1], [a_count a_count], '-ob', 'LineWidth', 5)
                xlim([0 100])
                ylim([0 42])

                hold on

                clear temp_fill
            end

        a_mega_matrix=[a_mega_matrix; a_vector];
        
        clear temp_burst_on temp_burst_dur temp_burst_off a_vector
        
    end
    
end

xlabel('Motif duration (%)')
ylabel('Cell ID')

a_summed=sum(a_mega_matrix);
a_summed_prob=smooth(a_summed./a_count); % smoothed using 5 point average

subplot(2,2,4)

plot(a_summed_prob, '-b', 'LineWidth', 2)
ylabel('Burst probability')
xlabel('Motif duration (%)')
xlim([0 100])
ylim([0 0.13])


[tbl, chi2, p]=crosstab(a_summed_prob, j_summed_prob)

clear tbl chi2 p

%% calculate mean single spikes per cell:

% adults:

for u=1:length(a_single_spikes)
   
    a_mean_single_spikes(u)=mean(a_single_spikes{1,u});
    
end

% juveniles:

for u=1:length(j_single_spikes)
   
    j_mean_single_spikes(u)=mean(j_single_spikes{1,u});
    
end

juv_scatter2=repmat(1,length(j_mean_single_spikes),1);

adt_scatter2=repmat(2,length(a_mean_single_spikes),1);

figure

scatter(juv_scatter2, j_mean_single_spikes, 'ro', 'Jitter', 'on')
hold on
scatter(adt_scatter2, a_mean_single_spikes, 'bo', 'Jitter', 'on')
hold on
plot([0.8 1.2], [mean(j_mean_single_spikes) mean(j_mean_single_spikes)], '-k')
hold on
plot([1.8 2.2], [mean(a_mean_single_spikes) mean(a_mean_single_spikes)], '-k')
xlim([0 3])
ylabel('Single Spikes per Trial')
axis square
box off

%% calculate single spike jitter:

% adults:

a_spike_jitter=[];

for u=1:length(a_single_spikes)
   
    temp_spikes=a_single_spikes{1,u};
    common_spike=mode(temp_spikes);
    
%    a_spike_jitter=[a_spike_jitter, abs(temp_spikes-common_spike)];
    a_spike_jitter{1,u}=abs(temp_spikes-common_spike);
    a_mean_spike_jitter(u)=mean(a_spike_jitter{1,u});
    
    
    clear temp_spikes common_spike
end

% juveniles:

j_spike_jitter=[];

for u=1:length(j_single_spikes)
   
    temp_spikes=j_single_spikes{1,u};
    common_spike=mode(temp_spikes);
    
%    j_spike_jitter=[j_spike_jitter, abs(temp_spikes-common_spike)];
    j_spike_jitter{1,u}=abs(temp_spikes-common_spike);
    j_mean_spike_jitter(u)=mean(j_spike_jitter{1,u});
    
    clear temp_spikes common_spike
end

juv_scatter3=repmat(1,length(j_mean_spike_jitter),1);

adt_scatter3=repmat(2,length(a_mean_spike_jitter),1);

figure

scatter(juv_scatter3, j_mean_spike_jitter, 'ro', 'Jitter', 'on')
hold on
scatter(adt_scatter3, a_mean_spike_jitter, 'bo', 'Jitter', 'on')
hold on
plot([0.8 1.2], [mean(j_mean_spike_jitter) mean(j_mean_spike_jitter)], '-k')
hold on
plot([1.8 2.2], [mean(a_mean_spike_jitter) mean(a_mean_spike_jitter)], '-k')
xlim([0 3])
ylabel('Single Spike Jitter')
axis square
box off

%% test correlations

for c=1:length(a_inst_fr_burst)
    
    a_mean_cell_inst_fr_burst(c)=nanmean(a_inst_fr_burst{1,c});
    
end

for c=1:length(j_inst_fr_burst)
    
    j_mean_cell_inst_fr_burst(c)=nanmean(j_inst_fr_burst{1,c});
    
end

% figure
% 
% plot(a_mean_cell_inst_fr_burst, a_cell_amp, 'bo')
% hold on
% plot(j_mean_cell_inst_fr_burst, j_cell_amp, 'ro')
% 
% st1=fitlm(j_mean_cell_inst_fr_burst, j_cell_amp)
% st2=fitlm(a_mean_cell_inst_fr_burst, a_cell_amp)
% 
% figure
% 
% plot(a_mean_cell_inst_fr_burst, a_freq_cells, 'bo')
% hold on
% plot(j_mean_cell_inst_fr_burst, j_freq_cells, 'ro')
% 
% figure
% 
% plot(a_cell_amp, a_freq_cells, 'bo')
% hold on
% plot(j_cell_amp, j_freq_cells, 'ro')


%% look at burst onset jitter:

% adults:

for t=1:length(a_cell_burst_time)

        temp_ms_on=a_cell_burst_time{1,t};
        temp_idx=a_burst_time_idx{1,t};
        str=unique(temp_idx);
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];

        if isempty(str)==0
        
                for z=1:length(str)

                    id=find(str(z)==temp_idx);
                    
                    temp_mode=mode(temp_ms_on(id));
                    temp_jitter=abs(temp_ms_on(id)-temp_mode);
                    burst_onset_jitter(z)=mean(temp_jitter);

                    clear id temp_mode temp_jitter

                end
                
            a_burst_mean_onset_jitter{1,t}=burst_onset_jitter;

        else
            continue
%             j_burst_mean_perc_on{1,t}=NaN;
%             j_burst_mean_perc_dur{1,t}=NaN;
            
        end
    clear temp_per_dur temp_per_on burst_perc_on_mean burst_perc_dur_mean temp_idx str burst_onset_jitter temp_ms_on
end

% juveniles:


for t=1:length(j_cell_burst_time)
    

        temp_ms_on=j_cell_burst_time{1,t};
        temp_idx=j_burst_time_idx{1,t};
        str=unique(temp_idx);
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];

        if isempty(str)==0
        
                for z=1:length(str)

                    id=find(str(z)==temp_idx);
                
                    
                    temp_mode=mode(temp_ms_on(id));
                    temp_jitter=abs(temp_ms_on(id)-temp_mode);
                    burst_onset_jitter(z)=mean(temp_jitter);

                    clear id temp_mode temp_jitter
                end
            j_burst_mean_onset_jitter{1,t}=burst_onset_jitter;
    
        else
            continue
%             j_burst_mean_perc_on{1,t}=NaN;
%             j_burst_mean_perc_dur{1,t}=NaN;
            
        end
    clear temp_per_dur temp_per_on burst_perc_on_mean burst_perc_dur_mean temp_idx str burst_onset_jitter temp_ms_on
end


a1=cell2mat(a_burst_mean_onset_jitter);
j1=cell2mat(j_burst_mean_onset_jitter);

figure

scatter(repmat(2,1,length(a1)),a1,'bo','Jitter', 'on')
hold on
scatter(repmat(1,1,length(j1)),j1,'ro','Jitter', 'on')
hold on
plot([0.8 1.2], [mean(j1) mean(j1)], '-k')
hold on
plot([1.8 2.2], [mean(a1) mean(a1)], '-k')
xlim([0 3])
axis square
box off
ylabel('Burst onset jitter (ms)')


%% plot within burst analysis:

figure

subplot(1,3,1)

scatter(repmat(1,1,length(j_burst_mean_mat)), j_burst_mean_mat, 'ro', 'Jitter', 'on')
hold on
scatter(repmat(2,1,length(a_burst_mean_mat)), a_burst_mean_mat, 'ko', 'Jitter', 'on')
hold on
plot([0.8 1.2], [mean(j_burst_mean_mat) mean(j_burst_mean_mat)], '-k')
hold on
plot([1.8 2.2], [mean(a_burst_mean_mat) mean(a_burst_mean_mat)], '-k')
xlim([0 3])
ylabel('Mean burst duration')
axis square


subplot(1,3,2)

scatter(repmat(1,1,length(j_burst_mean_mat)), j_mean_spikes_per_burst_mat, 'ro', 'Jitter', 'on')
hold on
scatter(repmat(2,1,length(a_burst_mean_mat)), a_mean_spikes_per_burst_mat, 'ko', 'Jitter', 'on')
hold on
plot([0.8 1.2], [mean(j_mean_spikes_per_burst_mat) mean(j_mean_spikes_per_burst_mat)], '-k')
hold on
plot([1.8 2.2], [mean(a_mean_spikes_per_burst_mat) mean(a_mean_spikes_per_burst_mat)], '-k')
xlim([0 3])
xlim([0 3])
ylim([1 10])
ylabel('Mean spikes per burst')
axis square

subplot(1,3,3)

% scatter(repmat(1,1,length(j_burst_mean_mat)), j_mean_burst_fr_mat, 'ro', 'Jitter', 'on')
% hold on
% scatter(repmat(2,1,length(a_burst_mean_mat)), a_mean_burst_fr_mat, 'ko', 'Jitter', 'on')
% hold on
% xlim([0 3])
% ylabel('Mean burst firing rate')
% axis square

j_burst_mean_isi=j_burst_mean_mat./(j_mean_spikes_per_burst_mat-1);

a_burst_mean_isi=a_burst_mean_mat./(a_mean_spikes_per_burst_mat-1);
% 
% figure

scatter(repmat(1,1,length(j_burst_mean_mat)), j_burst_mean_isi, 'ro', 'Jitter', 'on')
hold on
scatter(repmat(2,1,length(a_burst_mean_mat)), a_burst_mean_isi, 'ko', 'Jitter', 'on')
hold on
plot([0.8 1.2], [mean(j_burst_mean_isi) mean(j_burst_mean_isi)], '-k')
hold on
plot([1.8 2.2], [mean(a_burst_mean_isi) mean(a_burst_mean_isi)], '-k')
xlim([0 3])
ylabel('Mean burst ISI')
axis square

% look at burst duration per cell:

c1=0;

for j=1:length(a_burst_mean_duration)
   
    if isnan(a_burst_mean_duration{1,j})==0
    
        c1=c1+1;
        
    a_burst_dur(c1)=mean(a_burst_mean_duration{1,j});
    
    end
    
end

c1=0;

for j=1:length(j_burst_mean_duration)
   
    if isnan(j_burst_mean_duration{1,j})==0
    
        c1=c1+1;
        
    j_burst_dur(c1)=mean(j_burst_mean_duration{1,j});
    
    end
    
end


%% look at burst ISIs

% adults:


a_burst_isi_fr=[];
a_min_spikes_per_burst=[];

for t=1:length(a_burst_duration)
    
    figure

        temp_burst_dur=a_burst_duration{1,t};
        temp_isi=a_burst_isi{1,t};
        temp_idx=a_burst_time_idx_all{1,t};
        str=unique(a_burst_time_idx_all{1,t});
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];

        if isempty(str)==0
        
                for z=1:length(str)
                    
                    subplot(1,length(str),z)
                        axis square
                        box off

                    id=find(str(z)==temp_idx); % find where the burst idx is the same
                    
                    %for zirgs=1:length(id)
                            
                            temp_isi_mtrx=(temp_isi(id));
                            temp_id_burst_dur=temp_burst_dur(id);
                            
                            
                            for kaza=1:length(temp_isi_mtrx)
                               
                                nozeros=find(temp_isi_mtrx{1,kaza}>0);
                                
                                if isempty(nozeros)==0
                                
                                burst_isi_inst_fr(kaza,1:length(temp_isi_mtrx{1,kaza}(nozeros)))=1./temp_isi_mtrx{1,kaza}(nozeros)*1000; % calculates instantaneous fr for unique burst, expressed in sec
                                
                                plot(burst_isi_inst_fr(kaza,:), '-k')
                                hold on
                                
                                else
                                   
                                burst_isi_inst_fr(kaza,1:length(temp_isi_mtrx{1,kaza}))=1./temp_isi_mtrx{1,kaza}*1000;%./temp_id_burst_dur(kaza)*1000;    
                                
                                nozeros2=find(burst_isi_inst_fr(kaza,:)>0);
                                
                                plot(burst_isi_inst_fr(kaza,nozeros2), '-k')
                                hold on
                                
                                end
                                
                                clear nozeros nozeros2
                            end
                            
                            %temp_mat=cell2mat(burst_isi_inst_fr);
                            %temp_mean=mean(burst_isi_inst_fr);
                        
                            burst_mat_isi_fr{1,z}=burst_isi_inst_fr;
                            %min_spikes_per_burst(z)=min(a_spikes_per_burst{1,z}(str));
                            
                            clear burst_isi_inst_fr temp_isi_mtrx temp_id_burst_dur temp_mat temp_mean
                    
                    %end

                    clear id
                end

            a_burst_isi_fr=[a_burst_isi_fr, burst_mat_isi_fr];
           % a_min_spikes_per_burst=[a_min_spikes_per_burst, min_spikes_per_burst];
        else
            
            %a_burst_isi_fr{1,t}=NaN;
            
        end
    clear temp_burst_dur temp_isi min_spikes_per_burst temp_idx str burst_mat_isi_inst_fr burst_mat_isi_fr min_spikes_per_burst
    

    
    %pause
    close all
end

% juveniles:


j_burst_isi_fr=[];
j_min_spikes_per_burst=[];

for t=1:length(j_burst_duration)
    
    figure

        temp_burst_dur=j_burst_duration{1,t};
        temp_isi=j_burst_isi{1,t};
        temp_idx=j_burst_time_idx_all{1,t};
        str=unique(j_burst_time_idx_all{1,t});
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];

        if isempty(str)==0
        
                for z=1:length(str)
                    
                    subplot(1,length(str),z)

                    id=find(str(z)==temp_idx); % find where the burst idx is the same
                    
                    %for zirgs=1:length(id)
                            
                            temp_isi_mtrx=(temp_isi(id));
                            temp_id_burst_dur=temp_burst_dur(id);
                            
                            
                            for kaza=1:length(temp_isi_mtrx)
                               
                                nozeros=find(temp_isi_mtrx{1,kaza}>0);
                                
                                if isempty(nozeros)==0
                                
                                burst_isi_inst_fr(kaza,1:length(temp_isi_mtrx{1,kaza}(nozeros)))=1./temp_isi_mtrx{1,kaza}(nozeros)*1000; % calculates instantaneous fr for unique burst, expressed in sec
                                
                                plot(burst_isi_inst_fr(kaza,:), '-r')
                                hold on
                                axis square
                                box off
                                
                                else
                                   
                                burst_isi_inst_fr(kaza,1:length(temp_isi_mtrx{1,kaza}))=1./temp_isi_mtrx{1,kaza}*1000;%./temp_id_burst_dur(kaza)*1000;    
                                
                                nozeros2=find(burst_isi_inst_fr(kaza,:)>0);
                                
                                plot(burst_isi_inst_fr(kaza,nozeros2), '-r')
                                hold on
                                axis square
                                box off
                                
                                end
                                
                                clear nozeros nozeros2
                            end
                            
                            %temp_mat=cell2mat(burst_isi_inst_fr);
                            %temp_mean=mean(burst_isi_inst_fr);
                        
                            burst_mat_isi_fr{1,z}=burst_isi_inst_fr;
                           % min_spikes_per_burst(z)=min(j_spikes_per_burst{1,z}(str));
                            
                            clear burst_isi_inst_fr temp_isi_mtrx temp_id_burst_dur temp_mat temp_mean
                    
                    %end

                    clear id
                end

            j_burst_isi_fr=[j_burst_isi_fr, burst_mat_isi_fr];
          %  j_min_spikes_per_burst=[j_min_spikes_per_burst, min_spikes_per_burst];
    
        else
            
%             j_burst_isi_fr=[j_burst_isi_fr, Na];
            
        end
    clear temp_burst_dur temp_isi temp_idx str burst_mat_isi_inst_fr burst_mat_isi_fr min_spikes_per_burst
    

    
    %pause
    close all
end

%% plot adult and juvenile burst instantaneous rate progression:

adult_burst_matrix=zeros(67,10);

figure

subplot(1,2,1)

for t=1:length(a_burst_isi_fr)
   
    current_burst=a_burst_isi_fr{1,t};
%     
%     current_burst=mean(a_burst_isi_fr{1,t});
%     adult_burst_matrix(t,1:length(current_burst))=current_burst;
    
    for g=1:size(current_burst,1)
    
        nozeros=find(current_burst(g,:)>0);
        
        plot(current_burst(g,nozeros), 'o-k')
        hold on
        
    end
    
    clear current_burst nozeros
end

% mean_adb=mean(adult_burst_matrix);
% 
% plot(mean_adb, '-b', 'LineWidth', 2)
ylim([0 1000])
axis square
box off

juvenile_burst_matrix=zeros(9,10);

subplot(1,2,2)

for t=1:length(j_burst_isi_fr)
   
    current_burst=j_burst_isi_fr{1,t};
%      current_burst=mean(j_burst_isi_fr{1,t});
%      juvenile_burst_matrix(t,1:length(current_burst))=current_burst;
    
    for g=1:size(current_burst,1)
    
        nozeros=find(current_burst(g,:)>0);
        
        plot(current_burst(g,nozeros), 'o-r')
        hold on
        
    end
    
    clear current_burst nozeros
end

% mean_jdb=mean(juvenile_burst_matrix);
% 
% plot(mean_jdb, '-b', 'LineWidth', 2)
ylim([0 1000])
axis square
box off

%% plot burst instantaneous firing rate depending on how many mode spikes per burst

% first calculate mode for all spikes (not just repeated as before)


% adults:

a_burst_mode=[];

for t=1:length(a_spikes_per_burst)

        temp_spikes=a_spikes_per_burst{1,t};
        temp_idx=a_burst_time_idx_all{1,t};
        str=unique(a_burst_time_idx_all{1,t});
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];


        if isempty(str)==0
        
                for z=1:length(str)

                    id=find(str(z)==temp_idx);
                    
                    a_burst_mode=[a_burst_mode, mode(temp_spikes(id))];

                    clear id
                end           
        end
    clear temp_spikes temp_idx str u
end

% juveniles:

j_burst_mode=[];

for t=1:length(j_spikes_per_burst)

        temp_spikes=j_spikes_per_burst{1,t};
        temp_idx=j_burst_time_idx_all{1,t};
        str=unique(j_burst_time_idx_all{1,t});
        u=find(isnan(str)==1); % clear unique values of NaNs
        str(u)=[];


        if isempty(str)==0
        
                for z=1:length(str)

                    id=find(str(z)==temp_idx);
                    
                    j_burst_mode=[j_burst_mode, mode(temp_spikes(id))];

                    clear id
                end           
        end
    clear temp_spikes temp_idx str u
end


% now create an array of burst isi only for bursts with mode amount of
% spikes:



for g=1:length(a_burst_isi_fr)
    
    c1=0;

    temp_isi_fr=a_burst_isi_fr{1,g};
    temp_mode=a_burst_mode(g);
    
    
    for k=1:size(temp_isi_fr,1)
        
        idx=find(temp_isi_fr(k,:)>0);
        
        if length(idx)==(temp_mode-1) % isi will always be 1 less than spikes
          
            c1=c1+1;
            
            mode_burst_isi_fr(c1,:)=temp_isi_fr(k,idx);
   
        end
        
        clear idx
    end
    
    a_mode_burst_isi_fr{1,g}=mode_burst_isi_fr;
    a_mode_mean_burst_isi_fr{1,g}=mean(mode_burst_isi_fr);
    
    clear temp_mode temp_isi_fr mode_burst_isi_fr
end




for g=1:length(j_burst_isi_fr)
    
    c1=0;

    temp_isi_fr=j_burst_isi_fr{1,g};
    temp_mode=j_burst_mode(g);
    
    
    for k=1:size(temp_isi_fr,1)
        
        idx=find(temp_isi_fr(k,:)>0);
        
        if length(idx)==(temp_mode-1) % isi will always be 1 less than spikes
          
            c1=c1+1;
            
            mode_burst_isi_fr(c1,:)=temp_isi_fr(k,idx);
   
        end
        
        clear idx
    end
    
    j_mode_burst_isi_fr{1,g}=mode_burst_isi_fr;
    j_mode_mean_burst_isi_fr{1,g}=mean(mode_burst_isi_fr);
    
    clear temp_mode temp_isi_fr mode_burst_isi_fr
end


figure

subplot(1,2,1)

for t=1:length(a_mode_burst_isi_fr)
   
    plot(a_mode_mean_burst_isi_fr{1,t}, 'o-k')
    hold on
    
end
% 
ylim([0 1000])
ylabel('Instantaneous firing rate')
xlim([1 10])
xlabel('Number of spikes')
axis square
box off

subplot(1,2,2)

for t=1:length(j_mode_burst_isi_fr)
   
    plot(j_mode_mean_burst_isi_fr{1,t}, 'o-r')
    hold on
    
end

ylim([0 1000])
ylabel('Instantaneous firing rate')
xlim([1 10])
xlabel('Number of spikes')
axis square
box off

% plot only up to 5 spikes per burst for comparison:

figure

subplot(1,2,1)

a_5_spikes=[];
a_5_spikes_idx=[];

for t=1:length(a_mode_burst_isi_fr)
    
    if length(a_mode_mean_burst_isi_fr{1,t})<=5
   
        a_5_spikes=[a_5_spikes, a_mode_mean_burst_isi_fr{1,t}];
        a_5_spikes_idx=[a_5_spikes_idx, [1:length(a_mode_mean_burst_isi_fr{1,t})]];
        
        plot(a_mode_mean_burst_isi_fr{1,t}, 'o-k')
        hold on
    
    end
    
end
% 
% f3=fit(a_5_spikes_idx.',a_5_spikes.', 'poly2')
% plot(f3, '-b')

ylim([0 1000])
ylabel('Instantaneous firing rate')
xlim([0 6])
xlabel('Number of ISI')
axis square
box off

subplot(1,2,2)

j_5_spikes=[];
j_5_spikes_idx=[];

for t=1:length(j_mode_burst_isi_fr)
    
    if length(j_mode_mean_burst_isi_fr{1,t})<=5
        
        j_5_spikes=[j_5_spikes, j_mode_mean_burst_isi_fr{1,t}];
        j_5_spikes_idx=[j_5_spikes_idx, [1:length(j_mode_mean_burst_isi_fr{1,t})]];
        
   
        plot(j_mode_mean_burst_isi_fr{1,t}, 'o-r')
        hold on
    
    end
    
end

% f4=fit(j_5_spikes_idx.',j_5_spikes.', 'poly2')
% plot(f4, '-b')

ylim([0 1000])
ylabel('Instantaneous firing rate')
xlim([0 6])
xlabel('Number of ISI')
axis square
box off

%% see if mean burst ISI are corellated with mean psp duration per cell

for t=1:length(a_burst_isi)
   
    temp_isi=cell2mat(a_burst_isi{1,t});
    temp_idx=find(temp_isi>0);
    temp_idx2=find(a_psp_duration{1,t}>0)

    a_mean_cell_burst_isi(t)=mean(temp_isi(temp_idx));
    a_mean_psp_duration(t)=mean(a_psp_duration{1,t}(temp_idx2));
    
    clear temp_isi temp_idx temp_idx2 temp_fr
end

for t=1:length(j_burst_isi)
   
    temp_isi=cell2mat(j_burst_isi{1,t});
    temp_idx=find(temp_isi>0);
    temp_idx2=find(j_psp_duration{1,t}>0)

    
    j_mean_cell_burst_isi(t)=mean(temp_isi(temp_idx));
    j_mean_psp_duration(t)=mean(j_psp_duration{1,t}(temp_idx2));
    
    clear temp_isi temp_idx temp_idx2 temp_fr
end

figure

subplot(1,2,1)

plot(a_mean_cell_burst_isi,a_mean_psp_duration, 'ko')

subplot(1,2,2)

plot(j_mean_cell_burst_isi,j_mean_psp_duration, 'ro')

%% look at how many bursts are reoccuring in juveniles vs adults:

for t=1:length(a_burst_time_idx) % find reoccuring bursts
    
    temp_idx=find(isnan(a_burst_time_idx{1,t})==0);

a_bursts_per_cell(t)=length(unique(a_burst_time_idx{1,t}(temp_idx)));

clear temp_idx
end


for t=1:length(j_burst_time_idx) 
    
    temp_idx=find(isnan(j_burst_time_idx{1,t})==0);

j_bursts_per_cell(t)=length(unique(j_burst_time_idx{1,t}(temp_idx)));

clear temp_idx
end

% also look at all bursts in juveniles and adults:


for t=1:length(a_burst_time_idx_all) % find all bursts
    
    temp_idx=find(isnan(a_burst_time_idx_all{1,t})==0);

a_bursts_per_cell_all(t)=length(unique(a_burst_time_idx_all{1,t}(temp_idx)));

clear temp_idx
end


for t=1:length(j_burst_time_idx_all) 
    
    temp_idx=find(isnan(j_burst_time_idx_all{1,t})==0);

j_bursts_per_cell_all(t)=length(unique(j_burst_time_idx_all{1,t}(temp_idx)));

clear temp_idx
end


adult_percentage=sum(a_bursts_per_cell)/length(a_burst_mode);

juvenile_percentage=sum(j_bursts_per_cell)/length(j_burst_mode);

figure

bar(1,juvenile_percentage*100,'r')
hold on
bar(2,adult_percentage*100,'k')
ylim([0 100])
axis square
box off



%% correlate syllable ISI with neural dynamics:

close all



for z = 1:length(ISI_filelist)
    
    
name = char(ISI_filelist(z,1).name);



if z>10
    
    ISI_bird{1,z}=name(1:6);
    
elseif z>1 && z<5
    
    ISI_bird{1,z}=name(1:7);
    
elseif z==10
    
    ISI_bird{1,z}=name(1:6);
    
else
    
    ISI_bird{1,z}=name(1:4);
    
end

clear name 
end

% look at psp duration first:



for z=1:length(adu_bird)
   
   for g=1:length(ISI_bird)
    
       if strcmp(ISI_bird{1,g},adu_bird{1,z}) % see if bird names are the same
           
           a_mean_syll_ISI(z)=mean_syll_isi_ms(g); % in case names ar the same, save it as a similarity
           
           
       end
end 
end

for z=1:length(juv_bird)
   
   for g=1:length(ISI_bird)
    
       if strcmp(ISI_bird{1,g}(end-2:end),juv_bird{1,z}(end-2:end)) % see if bird names are the same
           
           j_mean_syll_ISI(z)=mean_syll_isi_ms(g); % in case names ar the same, save it as a similarity
           
           
       end
end 
end

% get rid of outliers:

j_drop=find(j_mean_psp_duration<60);
a_drop=find(a_mean_psp_duration<60);


figure

lr1=fitlm(a_mean_psp_duration(a_drop), a_mean_syll_ISI(a_drop))

lr2=fitlm(j_mean_psp_duration(j_drop), j_mean_syll_ISI(j_drop))

subplot(2,2,1)

plot(a_mean_psp_duration(a_drop), a_mean_syll_ISI(a_drop), 'ok')
axis square
box off
ylabel('Mean Syllable ISI (ms)')
xlabel('Mean PSP duration (ms)')
ylim([90 220])
xlim([10 55])

subplot(2,2,2)

plot(j_mean_psp_duration(j_drop), j_mean_syll_ISI(j_drop), 'ro')
axis square
box off
ylabel('Mean Syllable ISI (ms)')
xlabel('Mean PSP duration (ms)')
ylim([90 220])
xlim([10 55])

subplot(2,2,3)

plot(lr1)
ylabel('Mean Syllable ISI (ms)')
xlabel('Mean PSP duration (ms)')
ylim([90 220])
xlim([10 55])
axis square
box off

subplot(2,2,4)

plot(lr2)
ylabel('Mean Syllable ISI (ms)')
xlabel('Mean PSP duration (ms)')
ylim([90 220])
xlim([10 55])
axis square
box off


figure


lr3=fitlm(a_mean_syll_ISI, a_cell_amp)

lr4=fitlm(j_mean_syll_ISI, j_cell_amp)

subplot(2,2,1)

plot(a_mean_syll_ISI, a_cell_amp, 'ok')
axis square
box off
xlabel('Mean Syllable ISI (ms)')
ylabel('Mean PSP amplitude (mV)')
xlim([90 220])
ylim([2 9])

subplot(2,2,2)

plot(j_mean_syll_ISI, j_cell_amp, 'ro')
axis square
box off
xlabel('Mean Syllable ISI (ms)')
ylabel('Mean PSP amplitude (mV)')
xlim([90 220])
ylim([2 9])

subplot(2,2,3)

plot(lr3)
xlabel('Mean Syllable ISI (ms)')
ylabel('Mean PSP amplitude (mV)')
xlim([90 220])
ylim([2 9])
axis square
box off

subplot(2,2,4)

plot(lr4)
xlabel('Mean Syllable ISI (ms)')
ylabel('Mean PSP amplitude (mV)')
xlim([90 220])
ylim([2 9])
axis square
box off


figure

lr_psp_adult=fitlm(a_mean_psp_duration(a_drop), a_cell_amp(a_drop))

lr_psp_juvenile=fitlm(j_mean_psp_duration(j_drop), j_cell_amp(j_drop))

subplot(2,2,1)

plot(a_mean_psp_duration(a_drop),a_cell_amp(a_drop),'ko')
xlim([10 50])
ylim([2 9])
xlabel('Mean PSP duration (ms)')
ylabel('Mean PSP Amplitude (mV)')
axis square
box off

subplot(2,2,2)

plot(j_mean_psp_duration(j_drop),j_cell_amp(j_drop),'ro')
xlim([10 50])
ylim([2 9])
xlabel('Mean PSP duration (ms)')
ylabel('Mean PSP Amplitude (mV)')
axis square
box off

subplot(2,2,3)

plot(lr_psp_adult)
xlim([10 50])
ylim([2 9])
xlabel('Mean PSP duration (ms)')
ylabel('Mean PSP Amplitude (mV)')
axis square
box off

subplot(2,2,4)

plot(lr_psp_juvenile)
xlim([10 50])
ylim([2 9])
xlabel('Mean PSP duration (ms)')
ylabel('Mean PSP Amplitude (mV)')
axis square
box off



% figure burst ISI per cell

c1=0;

for t=1:length(a_bursts_per_cell_all)
   
    c1=c1+1;
    
    a_burst_mean_isi_per_cell{1,t}=a_burst_mean_isi(c1:c1+a_bursts_per_cell_all(t)-1);
    a_burst_mean_inst_fr_cell{1,t}=a_mean_burst_fr_mat(c1:c1+a_bursts_per_cell_all(t)-1);
    
    c1=sum(a_bursts_per_cell_all(1:t));
    
end


c1=0;

for t=1:length(j_bursts_per_cell_all)
   
    c1=c1+1;
    
    j_burst_mean_isi_per_cell{1,t}=j_burst_mean_isi(c1:c1+j_bursts_per_cell_all(t)-1);
    j_burst_mean_inst_fr_cell{1,t}=j_mean_burst_fr_mat(c1:c1+j_bursts_per_cell_all(t)-1);
    
    c1=sum(j_bursts_per_cell_all(1:t));
    
end

% plot the burst ISI correlation:

figure

subplot(2,2,1)

c1=0;
adu_burst_isi=[];
adu_syll_isi=[];
adu_burst_instfr=[];

for t=1:length(a_burst_mean_isi_per_cell)
   
    if isempty(a_burst_mean_isi_per_cell{1,t})==0 % make sure to only plot cells with bursts
    
        c1=c1+1;
        
%     plot(a_burst_mean_isi_per_cell{1,t}, a_mean_syll_ISI(t), 'ko')
%     hold on
    
    plot(a_burst_mean_inst_fr_cell{1,t}, a_mean_syll_ISI(t), 'ko')
    hold on
    
    adu_burst_instfr=[adu_burst_instfr, a_burst_mean_inst_fr_cell{1,t}];
    %adu_burst_isi=[adu_burst_isi,a_burst_mean_isi_per_cell{1,t}];
    adu_syll_isi=[adu_syll_isi, repmat(a_mean_syll_ISI(t),length(a_burst_mean_inst_fr_cell{1,t}),1).'];
    
    adu_mean_burst_isi(c1)=mean(a_burst_mean_isi_per_cell{1,t});
    adu_idx(c1)=t;    
    
    end
end
%plot(adu_burst_isi, adu_syll_isi, 'ko')
axis square
box off
ylabel('Mean Syllable ISI (ms)')
xlabel('Mean burst instantaneous firing rate (Hz)')
ylim([90 220])
xlim([0 800])

subplot(2,2,2)

c1=0;
juv_burst_isi=[];
juv_syll_isi=[];
juv_burst_instfr=[];

for t=1:length(j_burst_mean_isi_per_cell)
    
    if isempty(j_burst_mean_isi_per_cell{1,t})==0 
        
      c1=c1+1;
   
%     plot(j_burst_mean_isi_per_cell{1,t}, j_mean_syll_ISI(t), 'ro')
%     hold on
% 
    plot(j_burst_mean_inst_fr_cell{1,t}, j_mean_syll_ISI(t), 'ro')
    hold on
    
    
    juv_burst_instfr=[juv_burst_instfr, j_burst_mean_inst_fr_cell{1,t}];
    %juv_burst_isi=[juv_burst_isi, j_burst_mean_isi_per_cell{1,t}];
    juv_syll_isi=[juv_syll_isi, repmat(j_mean_syll_ISI(t),length(j_burst_mean_inst_fr_cell{1,t}),1).'];
    
    
    juv_mean_burst_isi(c1)=mean(j_burst_mean_isi_per_cell{1,t});
    juv_idx(c1)=t;
    
    end
end
%plot(juv_burst_isi, juv_syll_isi, 'ro')
axis square
box off
ylabel('Mean Syllable ISI (ms)')
xlabel('Mean burst instantaneous firing rate (Hz)')
ylim([90 220])
xlim([0 800])

lr5=fitlm(juv_burst_instfr, juv_syll_isi)

lr6=fitlm(adu_burst_instfr, adu_syll_isi)

subplot(2,2,3)

plot(lr6)
ylabel('Mean Syllable ISI (ms)')
xlabel('Mean burst instantaneous firing rate')
ylim([90 220])
xlim([0 800])
axis square
box off

subplot(2,2,4)

plot(lr5)
ylabel('Mean Syllable ISI (ms)')
xlabel('Mean burst instantaneous firing rate')
ylim([90 220])
xlim([0 800])
axis square
box off

%% exclude PSPs above 30 ms to avoid outliers:

for t=1:length(a_psp_duration)
   
    her=find(a_psp_duration{1,t}<30);
    
    a_psp_mean_duration_30ms(t)=mean(a_psp_duration{1,t}(her));
    
    clear her
end

for t=1:length(j_psp_duration)
   
    her=find(j_psp_duration{1,t}<30);
    
   j_psp_mean_duration_30ms(t)=mean(j_psp_duration{1,t}(her));
    
    clear her
end

figure

subplot(1,2,1)

plot(a_psp_mean_duration_30ms,a_cell_amp,'ko')
hold on

subplot(1,2,2)

plot(j_psp_mean_duration_30ms,j_cell_amp,'ro')


figure

plot(1,a_psp_mean_duration_30ms,'ko')
hold on
plot(2,j_psp_mean_duration_30ms,'ro')
xlim([0 3])

%% look at the corrected psp durations:

for t=1:length(a_psp_duration_corrected)
    
   a_psp_duration_corrected_mean(t)=nanmean(a_psp_duration_corrected{1,t});
   a_psp_amplitude_cells_corrected_mean(t)=nanmean(a_psp_amplitude_cells_corrected{1,t});
end

for t=1:length(j_psp_duration_corrected)
    
   j_psp_duration_corrected_mean(t)=nanmean(j_psp_duration_corrected{1,t});
   j_psp_amplitude_cells_corrected_mean(t)=nanmean(j_psp_amplitude_cells_corrected{1,t});
    
end

figure

plot(1,a_psp_duration_corrected_mean,'ko')
hold on
plot(2,j_psp_duration_corrected_mean,'ro')
xlim([0 3])

figure

subplot(1,2,1)

plot(j_psp_duration_corrected_mean,j_psp_amplitude_cells_corrected_mean,'ro')
axis square
ylim([2 9])
xlim([5 20])

subplot(1,2,2)

plot(a_psp_duration_corrected_mean,a_psp_amplitude_cells_corrected_mean,'ko')
axis square
ylim([2 9])
xlim([5 20])

%% look at number of spikes per trial mean:

figure

plot(1, j_nr_spikes_trial_mean, 'ro')
hold on
plot(2, a_nr_spikes_trial_mean, 'ko')
xlim([0 3])
axis square
box off

%% write models:

% create tables:

age=[repmat({str2},1,length(adu_bird)), repmat({str1},1,length(juv_bird))];

age_repeated_bursts=[repmat({str2},1,length(a_repeated_burst_bird)), repmat({str1},1,length(j_repeated_burst_bird))];

age_all_bursts=[repmat({str2},1,length(a_all_burst_bird)), repmat({str1},1,length(j_all_burst_bird))];

bird_table=table([adu_bird, juv_bird].', [a_cell_ID, j_cell_ID].', age.', [a_nr_spikes_trial_mean, j_nr_spikes_trial_mean].', [a_max_bursts, j_max_bursts].', [a_mean_single_spikes, j_mean_single_spikes].', 'VariableNames', {'bird_ID', 'cell_ID', 'age', 'nr_spikes_trial_mean', 'max_bursts', 'trial_mean_single_spikes'});

bird_repeated_burst_table=table([a_repeated_burst_bird, j_repeated_burst_bird].',[a_repeated_burst_cell_ID, j_repeated_burst_cell_ID].', age_repeated_bursts.', [a_jitters_mat, j_jitters_mat].', 'VariableNames', {'bird_ID', 'cell_ID', 'age', 'burst_jitter'});

bird_all_burst_table=table([a_all_burst_bird, j_all_burst_bird].', [a_all_burst_cell_ID, j_all_burst_cell_ID].',...
    age_all_bursts.', [a_burst_mean_isi, j_burst_mean_isi].', [a_burst_mean_mat, j_burst_mean_mat].',...
    [a_mean_spikes_per_burst_mat, j_mean_spikes_per_burst_mat].', [a_mean_burst_fr_mat, j_mean_burst_fr_mat].',...
    'VariableNames', {'bird_ID', 'cell_ID', 'age', 'burst_ISI', 'burst_duration', 'spikes_per_burst', 'burst_instfr'});

%


lme1=fitlme(bird_table,'nr_spikes_trial_mean~age+(1|bird_ID)+(1|bird_ID:cell_ID)')

[~,~,stats]=randomEffects(lme1);
lme2=fitlme(bird_table,'max_bursts~age+(1|bird_ID)+(1|bird_ID:cell_ID)')

[~,~,stats]=randomEffects(lme2);

lme3=fitlme(bird_table,'trial_mean_single_spikes~age+(1|bird_ID)+(1|bird_ID:cell_ID)')

[~,~,stats]=randomEffects(lme3);

lme4=fitlme(bird_repeated_burst_table,'burst_jitter~age+(1|bird_ID)+(1|bird_ID:cell_ID)')

[~,~,stats]=randomEffects(lme4);

lme5=fitlme(bird_all_burst_table,'burst_ISI~age+(1|bird_ID)+(1|bird_ID:cell_ID)')

lme6=fitlme(bird_all_burst_table,'burst_duration~age+(1|bird_ID)+(1|bird_ID:cell_ID)')

lme7=fitlme(bird_all_burst_table,'spikes_per_burst~age+(1|bird_ID)+(1|bird_ID:cell_ID)')

lme8=fitlme(bird_all_burst_table,'burst_instfr~age+(1|bird_ID)+(1|bird_ID:cell_ID)')