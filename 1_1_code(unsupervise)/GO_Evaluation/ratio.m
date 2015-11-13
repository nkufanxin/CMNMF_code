function [rat] = ratio(result)
%clc;
tic;
%% Load the matrix of sim(1274,1274)
load('similarity.mat');
%simi(simi==0) = 0.00001;
%simi = -10*log(simi);
% simi(simi==0) = 0.00001;
% simi=1./simi;
% simi = simi-diag(diag(simi));
%digits(4);
%sim = vpa(sim);
% simi = exp(simi/0.175);
% simi = simi-diag(diag(simi));
%% Load the matrix of cluster result(1274,1)
%load('cluster_result.mat');
%% Load the matrix of mgi id
load('mgi_id.mat');
% The cluster matrix of CMNMF
%result = result_CMNMF;
% The number of class
%total_cluster = size(result,2);
result(:,sum(result,1)==0)=[];
[m,n]=size(result)
total_cluster = size(result,2);
result = result';
[num_NonEmptyCluster,~] = size(result);
% Initial the matrix of inner class distance
inner_similarity = zeros(num_NonEmptyCluster,1);
%% Compute the distance of the inner class and save the result in indis
for i = 1:num_NonEmptyCluster
    %disp(['The average ',label,' inner distance of the ',num2str(i),'th class is computing']);
    % each row whose element is mgi id is one class and remove the zero element
    genes_in_each_cluster = result(i,result(i,:)>0);    
    % Initial the position matrix of mgi id
    %row = zeros(1,length(genes_in_each_cluster));
    % Put the position of mgi id in the row matrix   
    [~,row,~]= intersect(mgi_id,genes_in_each_cluster);    
    % the number of element in one class
    dim_row = length(row);
    % The pairs in the one class
    pairs = dim_row*(dim_row-1)/2;    
    % Compute the average distance of one class and save it in indis    
    dis = sum(sum(simi(row,row)))/2;   
    if pairs ~=0 % There are some classes whos element is one
        inner_similarity(i) = dis/pairs;
    end
end
%save(strcat(label,'.mat'),'indis');
%% Initial the matrix of the distance between two class
between_cluster_similarity = zeros(num_NonEmptyCluster,num_NonEmptyCluster);
% The sum of all betwdis
allbetw = 0;
for i = 1:num_NonEmptyCluster 
    % Find the first class and get the positon
    genes_in_each_cluster = result(i,result(i,:)>0);
    [~,ia,~]= intersect(mgi_id,genes_in_each_cluster);
    dis = 0;
    pairs_between_two_clusters = 0;
    % Find the second class and get the position
    for j = i+1:num_NonEmptyCluster
        %disp(['The average ',label,' distance between the ',num2str(i),'th class and the ',num2str(j),'th class is computing']);
        genes_in_each_cluster2 = result(j,result(j,:)>0);
        [~,ia2,~]= intersect(mgi_id,genes_in_each_cluster2); 
        pairs_between_two_clusters = length(ia)*length(ia2);
        if(pairs_between_two_clusters~=0)
        dis = sum(sum(simi(ia,ia2)));  
        between_cluster_similarity(i,j) = dis/pairs_between_two_clusters; 
        end
    end
end
%save(strcat(label,'.mat'),'indis','betwdis');
% The average distance of all classes in themselves
mean_inner_similarity = sum(inner_similarity)/total_cluster;
mean_between_similarity = sum(sum(between_cluster_similarity))/(total_cluster*(total_cluster-1)/2);
% The average distanve of all classes between each other

% The final ratio
rat = mean_inner_similarity/mean_between_similarity;
rat = 1/rat;
disp(['The Runtime of Program is:',num2str(toc)]);
end
