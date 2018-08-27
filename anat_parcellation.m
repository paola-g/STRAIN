addpath(genpath('NeuroElf_v09c'));
myvmp = xff('all_dmn.vmp');
load('subject_labels');
maps = myvmp.Map;
for i=1:140
    mynames{i}=maps(i).Name;
end
index = 0;
%%
medianvoi = zeros(115,90);
%meanvoi = zeros(115,90);
myvmp.NrOfMaps = 1;

for i=1:115
    index = strmatch(c{i},mynames);
    mymap = maps(index(1));
    myvmp.Map = mymap;
    myvtc = myvmp.ConvertToVTC();
    %bvoi = xff('Brodmann.voi');
    rvoi = xff('RVOI.voi');
    opts = struct;
    % one of 'cubic', 'lanczos3', 'linear', {'nearest'}
    opts.imeth = 'cubic';
    [medianvoi(i,:), voival, voin, wr,vois] = myvtc.VOITimeCourse(rvoi,opts);
end


% myvtc = xff('new:vtc');
% myvtc.NrOfVolumes = 1;
% data = myvmp.Map(i).VMPData;
% vtcdata = myvtc.VTCData;
% vtcdata(1,:,:,:) = data;
% myvtc.VTCData = vtcdata;

%% VOITimeCourseOrig
myvmp.NrOfMaps = 115;
for i=1:140
    mynames{i}=maps(i).Name;
end

for i=1:115
    index = strmatch(c{i},mynames);
    myindexes(i) = index(1);
    mymap = maps(index(1));
    newmaps(i) = mymap;
end

%bvoi = xff('Brodmann.voi');
rvoi = xff('RVOI.voi');
opts = struct;
% one of 'cubic', 'lanczos3', 'linear', {'nearest'}
opts.imeth = 'cubic';
myvmp.Map = newmaps;
%medianvois = zeros(115,90);

%vmp=xff('all_dmn.vmp') ;
%vtc=vmp.ConvertToVTC();
voi=xff('RVOI.voi');
vtc = myvmp.ConvertToVTC();
[medianvois,vois]=vtc.VOITimeCourseOrig(voi);

%%
load multi_disease_DMN_ext_fix DMNlabels
labels = DMNlabels;
control_idx = find(labels==1);
sla_idx = find(labels==2);
park_idx = find(labels==3);

%%
load medianvois

%% distribuzione voi per il primo soggetto

mymap = newmaps(1);
data = mymap.VMPData;
for i=1:90
hist(data(vois{i}),100);
mu=median(data(vois{i}));

%Overlay the median
hold on
plot([mu,mu],ylim,'r--','LineWidth',2)
hold off
pause
end

%% distribuzione voi per classe
data_all = zeros(115,58*40*46);
for i=1:115
   mymap = newmaps(i);
   data_all(i,:) = mymap.VMPData(:);
end


mean_data_all = mean(data_all);
mean_data_ctrl = mean(data_all(control_idx,:));
mean_data_sla = mean(data_all(sla_idx,:));
mean_data_park = mean(data_all(park_idx,:));

for j=1:19
i = uniontopk(j);
ax1 = subplot(2,2,1);
hist(mean_data_all(unique(vois{i})),100);
mu=median(mean_data_all(unique(vois{i})));
%Overlay the median
hold on
plot([mu,mu],ylim,'r--','LineWidth',2)
hold off
title('All')
ax2 = subplot(2,2,2);
hist(mean_data_ctrl(unique(vois{i})),100);
mu=median(mean_data_ctrl(unique(vois{i})));
%Overlay the median
hold on
plot([mu,mu],ylim,'r--','LineWidth',2)
hold off
title('Controls')
ax3 = subplot(2,2,3);
hist(mean_data_sla(unique(vois{i})),100);
mu=median(mean_data_sla(unique(vois{i})));
%Overlay the median
hold on
plot([mu,mu],ylim,'r--','LineWidth',2)
hold off
title('SLA')
ax4 = subplot(2,2,4);
hist(mean_data_park(unique(vois{i})),100);
mu=median(mean_data_park(unique(vois{i})));
%Overlay the median
hold on
plot([mu,mu],ylim,'r--','LineWidth',2)
hold off
title('Parkinson')

ylims = zeros(4,2);
ylims(1, :) = get(ax1,'ylim');
ylims(2, :) = get(ax2,'ylim');
ylims(3, :) = get(ax3,'ylim');
ylims(4, :) = get(ax4,'ylim');

xlims = zeros(4,2);
xlims(1, :) = get(ax1,'xlim');
xlims(2, :) = get(ax2,'xlim');
xlims(3, :) = get(ax3,'xlim');
xlims(4, :) = get(ax4,'xlim');

xmin = min(xlims(:,1));
ymin = min(ylims(:,1));
xmax = max(xlims(:,2));
ymax = max(ylims(:,2));
set([ax1 ax2 ax3 ax4],'xlim',[xmin xmax]);
set([ax1 ax2 ax3 ax4],'ylim',[ymin ymax]);
pause
saveas(gcf, sprintf('multidisease/median/histvoi_%d.png',uniontopk(j)));
saveas(gcf, sprintf('multidisease/median/histvoi_%d.fig',uniontopk(j)));
end

%%
load medianvois_control_topklists
topkcontrols = topkspace;
load medianvois_sla_topklists
topksla = topkspace;
load medianvois_park_topklists
topkpark = topkspace;
uniontopk = union(topkpark,topksla);
uniontopk = union(uniontopk,topkcontrols);

mean_medians_ctrl = mean(medianvois(control_idx, :));
medians_topkcontrols = mean_medians_ctrl(topkcontrols)';
mean_medians_sla = mean(medianvois(sla_idx, :));
medians_topksla = mean_medians_sla(topksla)';
mean_medians_park = mean(medianvois(park_idx, :));
medians_topkpark = mean_medians_park(topkpark)';

excludedmedians = mean_medians_ctrl;
excludedmedians(uniontopk) = [];
ax1 = subplot(1,2,1);
hist(medians_topkcontrols,10);
ax2 = subplot(1,2,2);
hist(excludedmedians,10);
minx = min(min([get(ax1, 'xlim') get(ax2, 'xlim')]));
maxx = max(max([get(ax1, 'xlim') get(ax2, 'xlim')]));
miny = min(min([get(ax2, 'ylim') get(ax2, 'ylim')]));
maxy = max(max([get(ax2, 'ylim') get(ax2, 'ylim')]));
set([ax1 ax2],'xlim',[minx maxx]);
set([ax1 ax2],'ylim',[miny maxy]);
%%
excluded_vois = vois;
excluded_vois(uniontopk) = [];

%%
for i=1:90
sizes(i) = length(unique(vois{i}));
end
sizes = sizes';
topsize = sizes(uniontopk);
hist(sizes,50)

%% input topklists
load medianvois
medians_ctrl = medianvois(control_idx, :);
medians_sla = medianvois(sla_idx, :);
medians_park = medianvois(park_idx, :);
[Y,rank_sla] = sort(medians_sla,2,'descend');
[Y,rank_park] = sort(medians_park,2,'descend');
[Y,rank_ctrl] = sort(medians_ctrl,2,'descend');
save md_anatvoi_desc_median_ranks rank_sla rank_park rank_ctrl
[Y,rank_sla] = sort(medians_sla,2,'ascend');
[Y,rank_park] = sort(medians_park,2,'ascend');
[Y,rank_ctrl] = sort(medians_ctrl,2,'ascend');
save md_anatvoi_asc_median_ranks rank_sla rank_park rank_ctrl
[Y,rank_sla] = sort(abs(medians_sla),2,'descend');
[Y,rank_park] = sort(abs(medians_park),2,'descend');
[Y,rank_ctrl] = sort(abs(medians_ctrl),2,'descend');
save md_anatvoi_abs_median_ranks rank_sla rank_park rank_ctrl

