function [label_final,N]=anatomicalClusterDistance(map,voxels,th,cluster_size)

% map Nx1
% voxels Nx3

pos_search_voxels=(find(abs(map)==th));
search_voxels=voxels(pos_search_voxels,:);

label_final=zeros(size(map));

if isempty(search_voxels)
    N=0;
    return
end

no_search=size(search_voxels,1);
label=zeros(no_search,1);
n=1;
clabel=0;
while n<=no_search;
    
    current_voxel=search_voxels(n,:);
    if label(n)==0
        clabel=clabel+1;
        label(n)=clabel;
        
        cand_x=find( (search_voxels(:,1)==current_voxel(1)-1) | (search_voxels(:,1)==current_voxel(1))  | (search_voxels(:,1)==current_voxel(1)+1));
        cand_y=find( (search_voxels(:,2)==current_voxel(2)-1) | (search_voxels(:,2)==current_voxel(2))  | (search_voxels(:,2)==current_voxel(2)+1));
        cand_z=find( (search_voxels(:,3)==current_voxel(3)-1) | (search_voxels(:,3)==current_voxel(3))  | (search_voxels(:,3)==current_voxel(3)+1));
        len=[length(cand_x) length(cand_y) length(cand_z)];
        m=min(find(len==min(len)));
        
        
        if m==1
            
            for i=1:length(cand_x)
                
                s=sum(cand_y==cand_x(i))+sum(cand_z==cand_x(i));
                if s==2
                    if label(cand_x(i))==0
                        label(cand_x(i))=clabel;
                    else   
                        cs=find(label==label(cand_x(i)));
                        label(cs)=clabel;
                    end %label
                end % s
                
            end  % for
        end  % if m
        
        
        if m==2
            
            for i=1:length(cand_y)
                
                s=sum(cand_x==cand_y(i))+sum(cand_z==cand_y(i));
                if s==2
                    if label(cand_y(i))==0
                        label(cand_y(i))=clabel;
                    else   
                        cs=find(label==label(cand_y(i)));
                        label(cs)=clabel;
                    end
                end
            end
        end  
        
        if m==3
            
            for i=1:length(cand_z)
                
                s=sum(cand_y==cand_z(i))+sum(cand_x==cand_z(i));
                if s==2
                    if label(cand_z(i))==0
                        label(cand_z(i))=clabel;
                    else   
                        cs=find(label==label(cand_z(i)));
                        label(cs)=clabel;
                    end
                end
            end
        end   
    end  % label
    
    n=n+1;
    
end  % while

%search_voxels

N=0;		 % number of active clusters
label_clu=zeros(size(label));
for i=1:max(label)
    
    nvox=find(label==i);
    
    if ~isempty(nvox)
        
        if length(nvox)>=cluster_size;
            
            label_clu(nvox)=ones(length(nvox),1);
            N=N+1;
        end
    end
end

label_final(pos_search_voxels)=label_clu;
