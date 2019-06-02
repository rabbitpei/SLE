clear;
clc;
close all;
%load alldata;
%load T1;
[T1] = load('data32_filled.txt');
load genename;
profile=T1;
%profile=zscore(profile');
%profile=profile';

fid=fopen('adj_edges_all.txt');
adjacent_network={};
j=0;
while ~feof(fid)
    tline=fgetl(fid);
    j=j+1;
    adjacent_network{j}=regexp(tline, '\t', 'split');
end
fclose(fid);
total_node_num=j;

%%
%profile=zscore(profile');
%profile=profile';
%profile(:,193)=mean(profile(:,192));
psize=size(profile);
for t=1:16
    for s=1:17
        tempcase(:,t,s)=profile(:,16*(s-1)+t);
    end
end
tempcontrol=tempcase(:,1:5,:);

psize=size(tempcase);
D=randperm(5);
%%
pretime=clock;
for l=1:psize(2)
    for na=1:total_node_num
        edges_list=[];
        center=adjacent_network{na}{1};
        e=0;
        for n=2:length(adjacent_network{na})
            nei=adjacent_network{na}{n};
            e=e+1;
            edges_list(e,:)=[str2num(center) str2num(nei)];
        end
        if e<1
            continue;
        end
        clear curr_pcc curr_pcc_add_onecase;
        
        for s=1:psize(3)
            for i=1:e
                curr_pcc(i)=abs(corr(reshape(tempcontrol(edges_list(i,1),D,s),5,1),...
                    reshape(tempcontrol(edges_list(i,2),D,s),5,1)));
                
                temp_add_onecase1=[reshape(tempcontrol(edges_list(i,1),D,s),1,5), ...
                    reshape(tempcase(edges_list(i,1),l,s),1,1)];
                temp_add_onecase2=[reshape(tempcontrol(edges_list(i,2),D,s),1,5),...
                    reshape(tempcase(edges_list(i,2),l,s),1,1)];
                curr_pcc_add_onecase(i)=abs(corr(temp_add_onecase1',temp_add_onecase2'));
            end
            
            prob_ctrl=curr_pcc/sum(curr_pcc);
            tmp_local_entropy=0;
            for i=1:e
                tmp_local_entropy=tmp_local_entropy-prob_ctrl(i)*log(prob_ctrl(i));
            end
            
            prob_add_onecase=curr_pcc_add_onecase/sum(curr_pcc_add_onecase);
            tmp_local_entropy_add_onecase=0;
            for i=1:e
                tmp_local_entropy_add_onecase=tmp_local_entropy_add_onecase-prob_add_onecase(i)*log(prob_add_onecase(i));
            end
            weighted_delta_entropy(na,l,s)=abs(tmp_local_entropy_add_onecase-tmp_local_entropy)*log(e);
        end
        if mod(na,1000)==0
            na,l
        end
    end
    
    currtime=clock;
    l,etime(currtime,pretime)
    pretime=currtime;
end
save weighted_delta_entropy weighted_delta_entropy;
load weighted_delta_entropy;
close all;
clear aver_weighted_delta_entropy_ctrl;
clear aver_weighted_delta_entropy_case;
ctrl_patients_idx=[2,3,4,9,11,14,16,17];
case_patients_idx=[1,5,6,7,8,10,12,13,15];
close all;
figure(1)
for s=1:length(ctrl_patients_idx)
    for l=1:psize(2)
        ss=ctrl_patients_idx(s);
        tmp_delta_entropy=sort(weighted_delta_entropy(:,l,ss),'descend');
        aver_weighted_delta_entropy_ctrl(l,s)=mean(tmp_delta_entropy(1:floor(total_node_num*1)));
    end
    plot([1:psize(2)],aver_weighted_delta_entropy_ctrl(1:psize(2),s),'b');
    hold on;
end
%title('Asymptomatic samples');
%figure(2)
for s=1:length(case_patients_idx)
    for l=1:psize(2)
        ss=case_patients_idx(s);
        tmp_delta_entropy=sort(weighted_delta_entropy(:,l,ss),'descend');
        aver_weighted_delta_entropy_case(l,s)=mean(tmp_delta_entropy(1:floor(total_node_num*1)));
    end
    %figure(s+1);
    hold on;
    plot([1:psize(2)],aver_weighted_delta_entropy_case(1:psize(2),s),'r');
    %title(['sample',num2str(ss)]);
    %hold off;
end
title('Symptomatic samples');

[se,idx]=sort(weighted_delta_entropy(:,7,2),'descend');
selected_deltaH_genes= symbols(idx(1:floor(total_node_num*1)));
fid=fopen('selected_deltaH_genes.txt','w');
for i=1:length(selected_deltaH_genes)
    fprintf(fid,'%s\n',selected_deltaH_genes{i});
end
fclose(fid);

