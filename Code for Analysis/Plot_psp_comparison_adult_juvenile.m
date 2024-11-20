clear all 
close all
[pathname] = uigetdir('DIRECTORY FOR FILES');
eval(['cd ' pathname]);

%% load data:

load('Adult_HVC_data_psp_duration_analysis2.mat')


a_cell_amp=average_amplitude_cells;

a_mean_cell_amp=average_amplitude_per_cell;

a_freq_cells=frequency_per_cell;

a_psp_mean_dur=psp_mean_duration_per_cell;


clearvars -except a_cell_amp a_freq_cells a_psp_mean_dur a_mean_cell_amp



load('Juvenile_HVC_data_psp_duration_analysis2.mat')

j_cell_amp=average_amplitude_cells;

j_mean_cell_amp=average_amplitude_per_cell;

j_freq_cells=frequency_per_cell;

j_psp_mean_dur=psp_mean_duration_per_cell;

clearvars -except a_cell_amp a_freq_cells j_cell_amp j_freq_cells a_psp_mean_dur j_psp_mean_dur a_mean_cell_amp...
    j_mean_cell_amp


%% plot figures:


max(a_cell_amp)
max(j_cell_amp)

edges=0:1:max(a_cell_amp)+1;

a_bin_idx=discretize(a_cell_amp,edges);
j_bin_idx=discretize(j_cell_amp,edges);

a_groups=unique(a_bin_idx);

for t=1:length(a_groups)
   
    a_howmany=find(a_groups(t)==a_bin_idx);
    a_percent(t)=(length(a_howmany)/length(a_cell_amp))*100;
    
    clear howmany
end

j_groups=unique(j_bin_idx);

for t=1:length(j_groups)
   
    j_howmany=find(j_groups(t)==j_bin_idx);
    j_percent(t)=(length(j_howmany)/length(j_cell_amp))*100;
    
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

edges2=[0,edges,max(edges)+0.5];

figure

stairs(edges2-0.5,a_mate2,'k')
hold on
stairs(edges2-0.5,j_mate2,'r')
hold on
plot([mean(j_cell_amp) mean(j_cell_amp)], [40 40], 'rv')
hold on
plot([mean(a_cell_amp) mean(a_cell_amp)], [40 40], 'kv')
ylabel('% PSPs')
xlabel('Amplitude (mV)')
xlim([0 14])
axis square
box off


[tbl, chi2, p]=crosstab(a_mate,j_mate);

clear tbl chi2 p

clear edges2 a_mate2 j_mate2 j_mate a_mate edges



figure

bar(a_groups,a_percent,'k')
hold on
bar(j_groups,j_percent,'r')
ylabel('% PSPs')
xlabel('Amplitude (mV)')
axis square
box off


figure
scatter(repmat(2,1,length(a_freq_cells)),a_freq_cells,'ko', 'Jitter', 'on')
hold on
scatter(repmat(1,1,length(j_freq_cells)),j_freq_cells,'ro', 'Jitter', 'on')
hold on
plot([1.8 2.2], [mean(a_freq_cells) mean(a_freq_cells)],'-k')
hold on
plot([0.8 1.2], [mean(j_freq_cells) mean(j_freq_cells)],'-r')
ytix2 = get(gca, 'YTick');
ylim([0 35])
xlim([0 3])
ylabel('Frequency (Hz)')
axis square
box off

%% plot psp duration, exclude outliers above 60 ms

j_drop=find(j_psp_mean_dur<60);
a_drop=find(a_psp_mean_dur<60);

figure

scatter(repmat(1,1,length(j_drop)), j_psp_mean_dur(j_drop), 'ro', 'Jitter', 'on')
hold on
scatter(repmat(2,1,length(a_drop)), a_psp_mean_dur(a_drop), 'ko', 'Jitter', 'on')
hold on
plot([0.8 1.2], [mean(j_psp_mean_dur(j_drop)) mean(j_psp_mean_dur(j_drop))],'-k')
hold on
plot([1.8 2.2], [mean(a_psp_mean_dur(a_drop)) mean(a_psp_mean_dur(a_drop))],'-k')
xlim([0 3])
ylabel('PSP duration (ms)')
axis square
box off

%% see if psp duration is correlated with amplitude:

figure

subplot(1,2,1)

plot(j_psp_mean_dur(j_drop), j_mean_cell_amp(j_drop), 'ro')
axis square

subplot(1,2,2)

plot(a_psp_mean_dur(a_drop), a_mean_cell_amp(a_drop), 'ko')
axis square

adult_lm=fitlm(a_psp_mean_dur(a_drop), a_mean_cell_amp(a_drop))

juv_lm=fitlm(j_psp_mean_dur(j_drop), j_mean_cell_amp(j_drop))