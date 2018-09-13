%% compute Jaccard index across subsamples
myfolder = '.';
matFiles = dir([myfolder,'consensus_toprank_d8-10-8v10_subranks*.mat']); % output of topklists_smp.R
nfiles = length(matFiles);
jdists_cl_smp = zeros(nfiles*(nfiles-1)/2, 3);
index = 1;
for j = 1:nfiles
for k = 1:nfiles
    if j<k
        try
        baseFileName_j = matFiles(j).name;
        fullFileName_j = fullfile(myfolder, baseFileName_j);
        load(fullFileName_j);
        ijfile = sscanf(baseFileName_j,'consensus_toprank_d8-10-8v10_subranks%d.mat'); % load output of md_kmeans_consensus_filter_subsamples
        load(sprintf('subsamples_filter_boot%d.mat',ijfile), 'cls_consensus', 'filter_idx')
        cl_filt_j = cls_consensus(filter_idx);
        cl_u_j = unique(cl_filt_j);
        cl_controls_j = cl_u_j(topkcontrols);
        cl_sla_j = cl_u_j(topksla);
        cl_park_j = cl_u_j(topkpark);
        vox_controls_j = filter_idx(find(ismember(cl_filt_j,cl_controls_j)));
        vox_sla_j = filter_idx(find(ismember(cl_filt_j,cl_sla_j)));
        vox_park_j= filter_idx(find(ismember(cl_filt_j,cl_park_j)));

        baseFileName_k = matFiles(k).name;
        fullFileName_k = fullfile(myfolder, baseFileName_k);
        load(fullFileName_k);
        ikfile = sscanf(baseFileName_k,'consensus_toprank_d8-10-8v10_subranks%d.mat');
        load(sprintf('subsamples_filter_boot%d.mat',ikfile), 'cls_consensus', 'filter_idx')
        cl_filt_k = cls_consensus(filter_idx);
        cl_u_k = unique(cl_filt_k);
        cl_controls_k = cl_u_k(topkcontrols);
        cl_sla_k = cl_u_k(topksla);
        cl_park_k = cl_u_k(topkpark);
        vox_controls_k = filter_idx(find(ismember(cl_filt_k,cl_controls_k)));
        vox_sla_k = filter_idx(find(ismember(cl_filt_k,cl_sla_k)));
        vox_park_k = filter_idx(find(ismember(cl_filt_k,cl_park_k)));

        % 1:controls 2:sla 3:park
        jdists_cl_smp(index, 1) = jaccard(vox_controls_j, vox_controls_k);
        jdists_cl_smp(index, 2) = jaccard(vox_sla_j, vox_sla_k);
        jdists_cl_smp(index, 3) = jaccard(vox_park_j, vox_park_k);
        index = index + 1;
        catch ME
            disp([j,k])
            continue
        end
    end
end
end

%% compute Jaccard index across subsamples (Jaccard between intersecting AAL vois)
load medianvois vois
topkspace_a = 1:90;

load mask4mm
maskIdx = find(mask4mm>0.9);

myfolder = '.';
matFiles = regexpdir(myfolder,'consensus_toprank_d8-10-8v10_subranks\d*.mat');
nfiles = length(matFiles);
jdists_cl_smp = zeros(nfiles*(nfiles-1)/2, 3);
index = 1;
for j = 1:nfiles
for k = 1:nfiles
    if j<k
        fullFileName_j = matFiles{j};
        tmp = splitstring(fullFileName_j,'\'); % change to '/' if unix
        baseFileName_j = tmp{end};
        load(fullFileName_j);
        ijfile = sscanf(baseFileName_j,'consensus_toprank_d8-10-8v10_subranks%d.mat');
        load(sprintf('subsamples_filter_boot%d.mat',ijfile), 'cls_consensus', 'filter_idx')
        cl_filt_j = cls_consensus(filter_idx);
        cl_u_j = unique(cl_filt_j);
        cl_controls_j = cl_u_j(topkcontrols);
        cl_sla_j = cl_u_j(topksla);
        cl_park_j = cl_u_j(topkpark);

        inters_controls_j = [];
        inters_sla_j = [];
        inters_park_j = [];
        for i=1:numel(topkcontrols)
            vol = zeros(size(mask4mm));
            a = vol(maskIdx(filter_idx));
            c = cl_filt_j==cl_controls_j(i);
            a(c) = 1;
            vol(maskIdx(filter_idx)) = a;
            volr = resampleVolume(vol, 2, 1/3,'nearest');
            idx_cl  = find(volr);
            for l=1:numel(topkspace_a)
                qnt = length(intersect(idx_cl, vois{topkspace_a(l)}));
                if(qnt > 0)
                   inters_controls_j = union(inters_controls_j,topkspace_a(l));
                end
            end    
        end
        for i=1:numel(topksla)
            vol = zeros(size(mask4mm));
            a = vol(maskIdx(filter_idx));
            c = cl_filt_j==cl_sla_j(i);
            a(c) = 1;
            vol(maskIdx(filter_idx)) = a;
            volr = resampleVolume(vol, 2, 1/3,'nearest');
            idx_cl  = find(volr);
            for l=1:numel(topkspace_a)
                qnt = length(intersect(idx_cl, vois{topkspace_a(l)}));
                if(qnt > 0)
                   inters_sla_j = union(inters_sla_j,topkspace_a(l));
                end
            end    
        end
        for i=1:numel(topkpark)
            vol = zeros(size(mask4mm));
            a = vol(maskIdx(filter_idx));
            c = cl_filt_j==cl_park_j(i);
            a(c) = 1;
            vol(maskIdx(filter_idx)) = a;
            volr = resampleVolume(vol, 2, 1/3,'nearest');
            idx_cl  = find(volr);
            for l=1:numel(topkspace_a)
                qnt = length(intersect(idx_cl, vois{topkspace_a(l)}));
                if(qnt > 0)
                   inters_park_j = union(inters_park_j,topkspace_a(l));
                end
            end    
        end


        
        
        fullFileName_k = matFiles{k};
        load(fullFileName_k);
        tmp = splitstring(fullFileName_k,'\');
        baseFileName_k = tmp{end};
        ikfile = sscanf(baseFileName_k,'consensus_toprank_d8-10-8v10_subranks%d.mat');
        load(sprintf('subsamples_filter_boot%d.mat',ikfile), 'cls_consensus', 'filter_idx')
        cl_filt_k = cls_consensus(filter_idx);
        cl_u_k = unique(cl_filt_k);
        cl_controls_k = cl_u_k(topkcontrols);
        cl_sla_k = cl_u_k(topksla);
        cl_park_k = cl_u_k(topkpark);

        inters_controls_k = [];
        inters_sla_k = [];
        inters_park_k = [];
        for i=1:numel(topkcontrols)
            vol = zeros(size(mask4mm));
            a = vol(maskIdx(filter_idx));
            c = cl_filt_k==cl_controls_k(i);
            a(c) = 1;
            vol(maskIdx(filter_idx)) = a;
            volr = resampleVolume(vol, 2, 1/3,'nearest');
            idx_cl  = find(volr);
            for l=1:numel(topkspace_a)
                qnt = length(intersect(idx_cl, vois{topkspace_a(l)}));
                if(qnt > 0)
                   inters_controls_k = union(inters_controls_k,topkspace_a(l));
                end
            end    
        end
        for i=1:numel(topksla)
            vol = zeros(size(mask4mm));
            a = vol(maskIdx(filter_idx));
            c = cl_filt_k==cl_sla_k(i);
            a(c) = 1;
            vol(maskIdx(filter_idx)) = a;
            volr = resampleVolume(vol, 2, 1/3,'nearest');
            idx_cl  = find(volr);
            for l=1:numel(topkspace_a)
                qnt = length(intersect(idx_cl, vois{topkspace_a(l)}));
                if(qnt > 0)
                   inters_sla_k = union(inters_sla_k,topkspace_a(l));
                end
            end    
        end
        for i=1:numel(topkpark)
            vol = zeros(size(mask4mm));
            a = vol(maskIdx(filter_idx));
            c = cl_filt_k==cl_park_k(i);
            a(c) = 1;
            vol(maskIdx(filter_idx)) = a;
            volr = resampleVolume(vol, 2, 1/3,'nearest');
            idx_cl  = find(volr);
            for l=1:numel(topkspace_a)
                qnt = length(intersect(idx_cl, vois{topkspace_a(l)}));
                if(qnt > 0)
                   inters_park_k = union(inters_park_k,topkspace_a(l));
                end
            end    
        end


        % 1:controls 2:sla 3:park
        jdists_cl_smp(index, 1) = jaccard(inters_controls_j, inters_controls_k);
        jdists_cl_smp(index, 2) = jaccard(inters_sla_j, inters_sla_k);
        jdists_cl_smp(index, 3) = jaccard(inters_park_j, inters_park_k);
        index = index + 1;
    end
end
end