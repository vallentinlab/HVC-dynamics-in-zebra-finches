%% If xcorr is NaN, check spike detection threshold for cutting spikes
clear all
close all
clc
%% Load data/ adress folder %% 
[pathname11] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname11]);

%% choose .mat file that includes aligned sound and recording traces

filelist2 = dir('*traces.mat');

number1 = length(filelist2);


%% variables:

fs = 40000; % Sampling rate

i= 1;

bird=[];
bird_ID = [];
treatment = [];
stdv_fr_call=[];
stdv_fr_playback=[];
bird_dph=[];
dph=[];
cell_group=[];

%% load data

for z = 1:number1
filename = char(strcat(pathname11,'\',filelist2(z,1).name));
load(filename)




%% for MD juvenile cells:

old_traces=traces;
stim_length=length(traces{1,1})/40;

number = size(traces,2);
interval=15; % interval in ms for looking at burst summation
thresh=20; % threshold for spike detection

trial_spikes_per_burst=[];
trial_burst_onsets=[];
cell_spikes_per_burst_max_var=[];
trial_burst_durations=[];
trial_burst_fr=[];
burst_isi=[];
trial_burst_pr_on=[];
trial_burst_pr_dur=[];
trial_single_spikes=[];
trial_burst_isi=[];
trial_burst_summation=[];
summation=[];


for k = 1:number

summation=[];
    
temp_avg=mean(traces{1,k}(:));
temp_trace=traces{1,k}(:)-temp_avg; % demean the traces here to cut spikes



%% Remove spikes (Jagadeesh 1997) %% removes spikes from all traces
spikes=[];
spikes= find(temp_trace>thresh);%{1,k}(:,1)>thresh);%-20); % gives you LOCS of where trace subthreshold is higher than thresh (spike)
for p = 1 : length(spikes)
        if spikes(p)<121 ||spikes(p)> size(traces{1,k}(:,1),1)-500 % || is an or condition, only executed if the part before is not true
        continue
    end
traces{1,k}([spikes(p)-120:spikes(p)+120],1)=NaN;
end

% see how many spikes per trial there are:

thresh_mVS=10;

if z==30 || z==54
    
    thresh_mVS=5;
    
elseif z==17
    
    thresh_mVS=15;
    
end

mVS=temp_trace-smooth(temp_trace,fs/80);

[pksp,time_spike_song] = findpeaks(mVS,'MinPeakHeight',thresh_mVS,'MinPeakDistance', 40); % changed to 40 (= 1 ms) for better detection
number_of_spikes_song(k) = length(findpeaks(mVS,'MinPeakHeight',thresh_mVS,'MinPeakDistance', 40));

%% plot a figure to make sure the threshold is correct:

% ef=figure('units','normalized','outerposition',[0 0 1 1]);
% 
% subplot(2,1,1)
% 
% plot(temp_trace,'-k')
% 
% subplot(2,1,2)
% 
% plot(mVS,'-r')
% hold on
% plot([0 length(temp_trace)], [thresh_mVS thresh_mVS], '-g')
% 
% pause
% close(ef)

%% see how many bursts you can find:

isi=(diff(time_spike_song))/40;

if z==17
    
    ids=find(isi<30); % adjust for the broad cell
    
elseif z==34

    ids=find(isi<3);
    
elseif z==18

    ids=find(isi<20);
    
elseif z==7

    ids=find(isi<10);
    
else

    ids=find(isi<15); % find where spikes are less than 20 ms apart, same interval as juveniles

end

bursts=0;
burst_time=[];

spikes_per_burst=0;
counter=1;

if length(ids)==1 % for a 2-spike burst
    
    bursts=1;
    spikes_per_burst(counter)=2;
    burst_time=(time_spike_song(ids(1)))/40; % express it in ms
    burst_off=(time_spike_song(ids(1)+1))/40;
    burst_duration(counter)=burst_off-burst_time(counter);
    burst_percent_on(counter)=burst_time(counter)/stim_length*100;
    burst_percent_dur(counter)=burst_duration(counter)/stim_length*100;
    burst_isi_dur{1,counter}=burst_duration; % this is in ms


elseif length(ids)>1 % define a burst
    
%     counter=counter+1;
    spikes_per_burst(counter)=2; % 1 burst means 2 spikes at least
    bursts=1; % if ids exceeds one, there is at least 1 burst
    burst_time=(time_spike_song((ids(1))))/40; % express it in ms
    burst_off=(time_spike_song(ids(1)+1))/40;
    burst_duration(counter)=burst_off-burst_time(counter);
    burst_percent_on(counter)=burst_time(counter)/stim_length*100;
    burst_percent_dur(counter)=burst_duration(counter)/stim_length*100;
    burst_isi_dur{1,counter}=burst_off-burst_time(counter);

    
    for o=2:length(ids)
   
        if ids(o)-ids(o-1)>1 % if not consecutive spikes, start a new burst

            bursts=bursts+1;
            counter=counter+1;
            spikes_per_burst(counter)=2;
            burst_time(counter)=(time_spike_song((ids(o))))/40;
            burst_off=(time_spike_song(ids(o)+1))/40; 
            burst_duration(counter)=burst_off-burst_time(counter);
            burst_percent_on(counter)=burst_time(counter)/stim_length*100;
            burst_percent_dur(counter)=burst_duration(counter)/stim_length*100;
            burst_isi_dur{1,counter}(o)=burst_off-burst_time(counter);


        
        elseif ids(o)-ids(o-1)==1 % if consecutive spikes, continue existing burst

            spikes_per_burst(counter)=spikes_per_burst(counter)+1;
            burst_off=(time_spike_song(ids(o)+1))/40; 
            burst_duration(counter)=burst_off-burst_time(counter);
            burst_percent_dur(counter)=burst_duration(counter)/stim_length*100;
            burst_isi_dur{1,counter}(o)=(time_spike_song((ids(o)+1))-time_spike_song((ids(o))))/40;
            
        end

    end

elseif length(ids)<1
    
    burst_isi_dur{1,counter}=NaN;
    burst_duration=0;
    bursts=0;
    burst_time=NaN;
    burst_percent_on=NaN;
    burst_percent_dur=NaN;
end

hp=find(spikes_per_burst>0); % exclude the no burst no spike events






trial_spikes_per_burst=[trial_spikes_per_burst, spikes_per_burst(hp)];
trial_bursts(k)=bursts;
trial_burst_onsets=[trial_burst_onsets, burst_time];
trial_burst_durations=[trial_burst_durations, burst_duration(hp)];
burst_inst_firing_rate=(spikes_per_burst(hp)-1)./(burst_duration(hp)/1000); % burst duration converted to s
trial_burst_fr=[trial_burst_fr, burst_inst_firing_rate];
trial_burst_pr_on=[trial_burst_pr_on, burst_percent_on];
trial_burst_pr_dur=[trial_burst_pr_dur, burst_percent_dur];
trial_burst_isi=[trial_burst_isi, burst_isi_dur(hp)];

trial_single_spikes=[trial_single_spikes, length(time_spike_song)-sum(spikes_per_burst)];

% gather together summation from each trial

figure

subplot(2,1,1)

plot(mVS)
hold on

for a=1:length(burst_time)
    
    if isnan(burst_time)==0
        
        if burst_time(a)>interval+1
        
        subplot(2,1,1)

        plot([burst_time(a)*40 (burst_time(a)+burst_duration(a))*40], [0.3 0.3], '-m', 'LineWidth',2)
        hold on
    
        summation{1,a}=temp_trace((burst_time(a)-(interval+1))*40:(burst_time(a)-1)*40); % leave out the 1 ms of half of the first spike
        
        subplot(2,1,2)

        plot(summation{1,a},'-r', 'LineWidth',1)
        hold on
        
        else
            
        summation{1,a}=NaN;
        
        end

    end
    
end

%pause
close all

trial_burst_summation=[trial_burst_summation, summation];

clear mVS pksp time_spike_song bursts isi ids burst_off summation



% interpolating traces

if isempty(spikes)==0 %check if the spikes were detected and cut
  traces_new{1,k}(:,1)=fixgaps(traces{1,k}(:,1));
  
else
    traces_new{1,k}(:,1)=traces{1,k}(:,1);
end  

% pause
% close(tr)


a=figure('units','normalized','outerposition',[0 0 1 1]);
% figure()

subplot(4,1,1)

plot(1/40:1/40:size(temp_trace,1)/40,temp_trace,'m')
hold on
xlim([1/40 size(temp_trace,1)/40])

mean_spikes(z)=mean(trial_spikes_per_burst);
std_spikes(z)=std(trial_spikes_per_burst);

for_text_bursts=num2str(trial_bursts(k));
for_text_spikes=num2str(spikes_per_burst(hp));
for_text_dur=num2str(burst_duration(hp));
for_text_fr=num2str(burst_inst_firing_rate);

% 
% text(size(temp_trace,1)/40/4,40,strcat('bursts per trial = ', for_text_bursts),'FontSize',12)
% 
 text(size(temp_trace,1)/40/2,40,strcat('spikes per burst = ', for_text_spikes),'FontSize',12)

%text(size(temp_trace,1)/40/2,40,strcat('burst duration = ', for_text_dur),'FontSize',12)

%text(size(temp_trace,1)/40/40,40,strcat('burst firing rate = ', for_text_fr),'FontSize',12)

ylim([-10 50])

subplot(4,1,2)

plot(1/40:1/40:size(traces_new{1,k},1)/40,traces_new{1,k},'k')
hold on
xlim([1/40 size(temp_trace,1)/40])


%pause
close all


clear temp_trace temp_avg spikes burst_inst_firing_rate burst_percent_on burst_percent_dur burst_isi_dur

clear for_text_bursts for_test_spikes_per_burst_mean for_test_spikes_per_burst_std hp spikes_per_burst burst_duration


end

% for_text_spikes_per_burst_mean=num2str(mean_spikes);
% for_text_spikes_per_burst_std=num2str(std_spikes);
% 
% subplot(4,1,3)
% 
% text(100,40,strcat('mean spikes per burst = ', for_text_spikes_per_burst_mean, '; std = ', for_text_spikes_per_burst_std),'FontSize',12)
% xlim([0 200])
% ylim([0 60])



clear mean_spikes_std_spikes
% clear spikes temp_trace temp_avg



%% continue the plot


% subplot(4,1,3)
% 
% cellmat=cell2mat(traces_new).';
% mean_vec=mean(cellmat);
% plot(1/40:1/40:length(mean_vec)/40,mean_vec,'b')
% xlim([1/40 length(mean_vec)/40])

%ram=figure(30)

%for t=1:number

% subplot(number*2,1,t)
% 
% spectrogram((sound{1,t}),1024,1000,1024,40000,'yaxis');
% colorbar off
% ylabel({'Frequency', '(Hz)'})% make it in two lines with cell array
% 
% myColorMap = jet; % Make a copy of jet.
% 
% % Assign black(all 0's) to black (the first row in myColorMap).
% myColorMap(1, :) = [0 0 0];
% colormap(myColorMap)
% ylim([1 7])
% caxis([-50 0])
% sound_envelope=smooth(abs(hilbert(sound{1,t})),50);

% clear cellmat 
% 
% subplot(number*2,1,number+t)
% 
% plot(1/40:1/40:size(traces_new{1,t},1)/40,traces_new{1,t},'k')
% hold on
% xlim([1/40 length(mean_vec)/40])
% xlim([1/40 size(traces_new{1,t},1)/40])

%end



%% De-mean data by substracting the mean value of each trace %% demeans all traces
for r = 1:number
averagemembranepotential_song = mean(traces_new{1,r},1);

demean_data_song(r,:) = traces_new{1,r}-averagemembranepotential_song;

traces_new_song{1,r}=demean_data_song(r,:).';


end


clear k


%% Calculate xcorr

traces_new_song_mat=cell2mat(traces_new_song).';

coefficient_son=[];
coefficient_song=[];

for m = 1:number
    
    for n = 1:number
        
        if n == m || n<m
            continue
        end
        %call:
        crosscorrall_song = xcorr(traces_new_song_mat(m,:),traces_new_song_mat(n,:),0,'coeff'); %40 000 stands for silent period (1s) before and after stimulus
        correlationcoefficients = crosscorrall_song;
        coefficient_son = [coefficient_son correlationcoefficients];

    end
end

% call:



% mean_song_trace{z,1}=mean_vec;

cell_burst_summation{1,z}=trial_burst_summation;

coefficient_song=  coefficient_son;
avg_corr_song = mean(coefficient_song);
std_corr_song = std(coefficient_song);
median_corr_song(z)=median(coefficient_song);

coefficient_song_all{z,1}=coefficient_song;
avg_corr_song_all(z)=avg_corr_song;
std_corr_song_all(z)=std_corr_song;

nr_bursts_trial{1,z}=trial_bursts;
nr_bursts_trial_mean(z)=mean(trial_bursts);
nr_bursts_trial_std(z)=std(trial_bursts);
nr_bursts_max(z)=max(trial_bursts);
nr_bursts_median(z)=median(trial_bursts);

nr_spikes_trial_cell{1,z}=number_of_spikes_song;
nr_spikes_trial_median(z)=median(number_of_spikes_song);
nr_spikes_trial_mean(z)=mean(number_of_spikes_song);
nr_spikes_trial_std(z)=std(number_of_spikes_song);
cell_group=[cell_group;repmat(z,length(number_of_spikes_song),1)];

cell_single_spikes{1,z}=trial_single_spikes;

if isempty(trial_spikes_per_burst)==1
   trial_spikes_per_burst=NaN;
   cell_spikes_per_burst{1,z}=NaN;
   cell_spikes_per_burst_variance(z)=NaN;
   cell_spikes_per_burst_max_var(z)=NaN;
   cell_spikes_per_burst_std(z)=NaN;
   
elseif length(trial_spikes_per_burst)==1

cell_burst_isi_duration{1,z}=trial_burst_isi;


    
cell_spikes_per_burst{1,z}=trial_spikes_per_burst;
cell_spikes_per_burst_variance(z)=0;
cell_spikes_per_burst_max_var(z)=NaN;
cell_spikes_per_burst_std(z)=std(trial_spikes_per_burst);

skipper=find(isnan(trial_burst_onsets)==0);

cell_burst_percent_onset{1,z}=trial_burst_pr_on(skipper);

cell_burst_percent_duration{1,z}=trial_burst_pr_dur(skipper);

cell_burst_time{1,z}=trial_burst_onsets(skipper);
cell_burst_duration{1,z}=trial_burst_durations;
cell_burst_firing_rate{1,z}=trial_burst_fr;

cell_burst_inst_firing_rate{1,z}=trial_burst_fr;
   
elseif length(trial_spikes_per_burst)>1
    
cell_burst_isi_duration{1,z}=trial_burst_isi;

cell_spikes_per_burst{1,z}=trial_spikes_per_burst;
cell_spikes_per_burst_variance(z)=var(trial_spikes_per_burst);
cell_spikes_per_burst_max_var(z)=max(diff(trial_spikes_per_burst));
cell_spikes_per_burst_std(z)=std(trial_spikes_per_burst);

skipper=find(isnan(trial_burst_onsets)==0);

cell_burst_percent_onset{1,z}=trial_burst_pr_on(skipper);

cell_burst_percent_duration{1,z}=trial_burst_pr_dur(skipper);

cell_burst_time{1,z}=trial_burst_onsets(skipper);
cell_burst_duration{1,z}=trial_burst_durations;
cell_burst_inst_firing_rate{1,z}=trial_burst_fr;

end
clearvars -except coefficient_call_all avg_corr_call_all std_corr_call all coefficient_play_all avg_corr_play_all...
    std_corr_play_all coefficient_silence_all avg_corr_silence_all std_corr_silence_all z pathname11 filelist2...
    fields number1 fs pre_silence post_silence sil bird bird_ID treatment treat dp dph cell cell_ID...
    pretime posttime i library thresh1 mean_call_trace dph bird_dph coefficient_song_all avg_corr_song_all...
    std_corr_song_all nr_spikes_trial_cell nr_spikes_trial_mean nr_spikes_trial_std nr_bursts_trial nr_bursts_trial_mean...
    nr_bursts_trial_std nr_bursts_max cell_group cell_spikes_per_burst cell_spikes_per_burst_variance cell_spikes_per_burst_std...
    cell_burst_time cell_spikes_per_burst_max_var mean_spikes std_spikes cell_burst_duration cell_burst_inst_firing_rate cell_burst_isi...
    cell_burst_percent_onset cell_burst_percent_duration cell_single_spikes cell_burst_isi_duration median_corr_song avg_corr_song...
    nr_spikes_trial_median nr_bursts_median cell_burst_summation

end


%% find corresponding bursts:
% burst is considered the same within 20 ms boundary

% for g=1:length(cell_burst_time)
%    
%     who=(cell_burst_time{1,g});
%     
%     for h=2:length(cell_burst_time{1,g})
%        
%         bet(h-1)=cell_burst_time{1,g}(h)-cell_burst_time{1,g}(h-1);
%         
%         if bet(h-1)<20
%             
%             c=c+1;
%             
%             idx(c)=h;
%         end
%         
%     end
%     
% end
% 
% a1=[1,3,5,7];
% a2=[1,4];
% a3=[2,5];
% a4=[2,4];
% a5=[3,5];
% 
% bursts=[{cell_spikes_per_burst{1,3}}, {cell_spikes_per_burst{1,4}}, {cell_spikes_per_burst{1,5}(a1)}, {cell_spikes_per_burst{1,6}}, {cell_spikes_per_burst{1,8}(a2)}, {cell_spikes_per_burst{1,8}(a3)}, {cell_spikes_per_burst{1,9}(a4)}, {cell_spikes_per_burst{1,9}(a5)}];
% 
% for h=1:length(bursts)
%    
%     burst_variance(h)=var(bursts{1,h});
%     
% end
% 
% 
% 
% %% plot figures:
% 
% sm=unique(burst_variance);
% 
% for t=1:length(sm)
%    
%     howmany=find(sm(t)==burst_variance);
%     percent(t)=(length(howmany)/length(burst_variance))*100;
%     
%     clear howmany
% end
% 
% figure
% 
% bar(sm,percent,'k')
% ylabel('% bursts')
% xlabel('burst variance')
% axis square
% box off
% 
% % or use this one:
% 
% histogram(burst_variance,length(sm))
% ytix = get(gca, 'YTick')
% set(gca, 'YTick',ytix, 'YTickLabel',ytix/length(burst_variance)*100) % change y axis to percentage here
% 
% clear percent
% 
% figure
% a=repmat(1,length(avg_corr_song_all),1);
% 
% scatter(a,avg_corr_song_all,'mo','Jitter','on')
% hold on
% plot([0.8 1.2], [mean(avg_corr_song_all) mean(avg_corr_song_all)], '-k')
% ylim([0 1])
% xlim([0 2])
% axis square
% box off
% ylabel('xcorr BOS')
% 
% u=[1:length(avg_corr_song_all)];
% 
% figure
% 
% plot(u, nr_spikes_trial_mean,'go')
% hold on
% errorbar(u, nr_spikes_trial_mean, nr_spikes_trial_std)
% hold on
% plot([0 11], [mean(nr_spikes_trial_mean) mean(nr_spikes_trial_mean)], '--k')
% xlim([0 11])
% axis square
% ylabel('Mean Nr of Spikes per trial')
% xlabel('Cell ID')
% 
% %plot nr_spikes boxplot:
% 
% spikes_trial=cell2mat(nr_spikes_trial_cell).';
% 
% cell_group2=cell_group*2;
% 
% figure
% 
% boxplot(spikes_trial,cell_group2)
% hold on
% beeswarm(cell_group,spikes_trial,'MarkerFaceColor','g')
% xlabel('Cell ID')
% ylabel('# spikes')
% axis square
% box off
% 
% figure
% 
% for z=1:length(nr_spikes_trial_cell)
% 
% scatter(repmat(z,length(nr_spikes_trial_cell{1,z}),1), nr_spikes_trial_cell{1,z},'bo','Jitter','on')
% hold on
% 
% end
% %errorbar(u, nr_bursts_trial_mean, nr_spikes_trial_std)
% %hold on
% plot([0 10], [mean(nr_spikes_trial_mean) mean(nr_spikes_trial_mean)], '--k')
% xlim([0 10])
% %ylim([0 4])
% axis square
% ylabel('Nr of Spikes per trial')
% xlabel('Cell ID')
% 
% figure
% 
% for z=1:length(nr_bursts_trial)
% 
% scatter(repmat(z,length(nr_bursts_trial{1,z}),1), nr_bursts_trial{1,z},'bo', 'Jitter','on')
% hold on
% 
% end
% %errorbar(u, nr_bursts_trial_mean, nr_spikes_trial_std)
% %hold on
% plot([0 10], [mean(nr_bursts_trial_mean) mean(nr_bursts_trial_mean)], '--k')
% xlim([0 10])
% ylim([0 4])
% axis square
% ylabel('Mean Nr of Bursts per trial')
% xlabel('Cell ID')
% 
% figure
% 
% plot(u, nr_bursts_trial_mean,'bo')
% hold on
% errorbar(u, nr_bursts_trial_mean, nr_bursts_trial_std)
% %hold on
% plot([0 11], [mean(nr_bursts_trial_mean) mean(nr_bursts_trial_mean)], '--k')
% xlim([0 11])
% axis square
% ylabel('Mean Nr of Bursts per trial')
% xlabel('Cell ID')
% 
% % plot the max nr of bursts:
% 
% a=unique(nr_bursts_max);
% 
% for t=1:length(a)
%    
%     howmany=find(a(t)==nr_bursts_max);
%     percent(t)=(length(howmany)/length(nr_bursts_max))*100;
%     
%     clear howmany
% end
% 
% figure
% 
% bar(a,percent,'r')
% ylim([0 100])
% ylabel('% cells')
% xlabel('# max bursts')
% axis square
% box off
% 
% figure
% bar(u,cell_spikes_per_burst_std, 'g')
% ylabel('std spikes per burst')
% xlabel('Cell_ID')
% axis square
% box off
% 
% figure
% bar(u,cell_spikes_per_burst_variance, 'b')
% ylabel('variance spikes per burst')
% xlabel('Cell_ID')
% axis square
% box off
% 
% 
% 
% trip=hist(cell_spikes_per_burst_variance);
% trip2=trip/length(cell_spikes_per_burst)*100;
% 
% figure
% bar(trip,trip2,'k')
% 
% burst_group=[];
% 
% for t=1:length(bursts)
% 
% burst_group=[burst_group; repmat(t,length(bursts{1,t}),1)];
% 
% end
% 
% burst_mat=cell2mat(bursts).';
% 
% figure
% 
% boxplot(burst_mat,burst_group)
% hold on
% beeswarm(burst_group,burst_mat)
% xlabel('Burst ID')
% ylabel('# spikes per burst')
% axis square
% box off