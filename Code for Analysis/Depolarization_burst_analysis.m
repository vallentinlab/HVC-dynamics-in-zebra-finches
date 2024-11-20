clear all 
close all
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);

%%
% take 15 ms before burst for subthresh analysis
filelist = dir('*.mat');

fs=40000;
thresh_mVS=0.1;

trial_spikes_per_burst=[];
trial_burst_onsets=[];
trial_burst_durations=[];
trial_burst_fr=[];
trial_burst_isi=[];
trial_single_spikes=[];
interval=15;

for k=1:length(filelist)

filename = char(strcat(pathname,'\',filelist(k,1).name));
load(filename)


trial_spikes_per_burst=[];
trial_burst_onsets=[];
trial_burst_durations=[];
trial_burst_fr=[];
trial_burst_isi=[];
trial_single_spikes=[];
interval=15;

%demean trace:

data_mean=mean(data_new);
data=data_new-data_mean;


mVS=data-smooth(data,fs/80);

[pksp,time_spike_song] = findpeaks(mVS,'MinPeakHeight',thresh_mVS,'MinPeakDistance', 40); % changed to 40 (= 1 ms) for better detection
number_of_spikes_song(k) = length(findpeaks(mVS,'MinPeakHeight',thresh_mVS,'MinPeakDistance', 40));

%% see how many bursts you can find:

isi=(diff(time_spike_song))/40;

ids=find(isi<100); % find where spikes are less than 100 ms apart, since you only work with depolarized events

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
    burst_isi_dur{1,counter}=burst_off-burst_time(counter);

    
    for o=2:length(ids)
   
        if ids(o)-ids(o-1)>1 % if not consecutive spikes, start a new burst

            bursts=bursts+1;
            counter=counter+1;
            spikes_per_burst(counter)=2;
            burst_time(counter)=(time_spike_song((ids(o))))/40;
            burst_off=(time_spike_song(ids(o)+1))/40; 
            burst_duration(counter)=burst_off-burst_time(counter);
            burst_isi_dur{1,counter}(o)=burst_off-burst_time(counter);


        
        elseif ids(o)-ids(o-1)==1 % if consecutive spikes, continue existing burst

            spikes_per_burst(counter)=spikes_per_burst(counter)+1;
            burst_off=(time_spike_song(ids(o)+1))/40; 
            burst_duration(counter)=burst_off-burst_time(counter);
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
trial_burst_isi=[trial_burst_isi, burst_isi_dur(hp)];

trial_single_spikes=[trial_single_spikes, length(time_spike_song)-sum(spikes_per_burst)];

figure

subplot(2,1,1)

plot(mVS)
hold on
for a=1:length(trial_burst_onsets)
    
subplot(2,1,1)

plot([trial_burst_onsets(a)*40 (trial_burst_onsets(a)+trial_burst_durations(a))*40], [0.3 0.3], '-m', 'LineWidth',2)
hold on

summation{1,a}=data((trial_burst_onsets(a)-interval)*40:trial_burst_onsets(a)*40);

subplot(2,1,2)

plot(summation{1,a},'-r', 'LineWidth',1)
hold on
end
% plot([1 length(mVS)], [thresh_mVS thresh_mVS], '-r', 'LineWidt', 2) 
% 
pause
close all

number=length(summation);

summation_mat=cell2mat(summation).';

coefficient_sum=[];
%coefficient_song=[];

for m = 1:number
    
    for n = 1:number
        
        if n == m || n<m
            continue
        end

        crosscorrall_summation = xcorr(summation_mat(m,:),summation_mat(n,:),0,'coeff');
        correlationcoefficients = crosscorrall_summation;
        coefficient_sum = [coefficient_sum correlationcoefficients];
        
    end
end

avg_corr_summation(k) = mean(coefficient_sum);
std_corr_summation(k) = std(coefficient_sum);

clearvars -except avg_corr_summation std_corr_summation filelist interval pathname fs k thresh_mVS

end