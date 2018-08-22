clc;
clear;
close all;

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

fpi=fopen('pruning_normal_expression_step2.txt');
hline = textscan(fpi, '%s', 1, 'delimiter', '\n');
field=textscan(hline{1}{1},'%s');
clear format;
format='%s';
% format=[format,' %s'];
for i=2:59
    format=[format,' %f'];
end
lines =textscan(fpi, format,1000000,'delimiter', '\t');
pipi=lines{1};
pprofile = [];
for i = 2 :59
    pprofile = [pprofile, lines{i}];
end
fclose(fpi);

fpi=fopen('pruning_tumor_expression_step2.txt');
hline_patient_IDs = textscan(fpi, '%s', 1, 'delimiter', '\n');
hline = textscan(fpi, '%s', 1, 'delimiter', '\n');
field=textscan(hline{1}{1},'%s');
clear format;
format='%s';
% format=[format,' %s'];
for i=2:426
    format=[format,' %f'];
end
lines =textscan(fpi, format,1000000,'delimiter', '\t');
mipi=lines{1};
mprofile = [];
for i = 2 :426
    mprofile = [mprofile, lines{i}];
end
fclose(fpi);

psize=size(pprofile);
allprofile(:,1:size(pprofile,2))=pprofile;
allprofile(:,size(pprofile,2)+1:size(pprofile,2)+size(mprofile,2))=mprofile;
%allprofile=zscore(allprofile');
%allprofile=allprofile';
pprofile=allprofile(:,1:size(pprofile,2));
mprofile=allprofile(:,size(pprofile,2)+1:size(pprofile,2)+size(mprofile,2));
tempcontrol=pprofile;
tempcase=zeros(psize(1),8,125);
patients_num=[3,106,124,39,59,62,10,21];
tempcase(:,1,1:patients_num(1))=mprofile(:,1:3);    % Stage I
tempcase(:,2,1:patients_num(2))=mprofile(:,4:109);   % Stage IA
tempcase(:,3,1:patients_num(3))=mprofile(:,110:233);  % Stage IB
tempcase(:,4,1:patients_num(4))=mprofile(:,234:272);   % Stage IIA
tempcase(:,5,1:patients_num(5))=mprofile(:,273:331);    % Stage IIB
tempcase(:,6,1:patients_num(6))=mprofile(:,332:393);   % Stage IIIA
tempcase(:,7,1:patients_num(7))=mprofile(:,394:403);    % Stage IIIB
tempcase(:,8,1:patients_num(8))=mprofile(:,404:424);     %IV
psize=size(tempcase);

%%
pretime=clock;
weighted_delta_entropy=zeros(total_node_num,psize(2),psize(3));
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
        for i=1:e
            curr_pcc(i)=abs(corr(tempcontrol(edges_list(i,1),:)',...
                tempcontrol(edges_list(i,2),:)'));
        end
        prob_ctrl=curr_pcc/sum(curr_pcc);
        tmp_local_entropy=0;
        for i=1:e
            tmp_local_entropy=tmp_local_entropy-prob_ctrl(i)*log(prob_ctrl(i));
        end
        
        for s=1:patients_num(l)    % psize(3)
            for i=1:e
                temp_add_onecase1=[tempcontrol(edges_list(i,1),:),...
                    reshape(tempcase(edges_list(i,1),l,s),1,1)];
                temp_add_onecase2=[tempcontrol(edges_list(i,2),:),...
                    reshape(tempcase(edges_list(i,2),l,s),1,1)];
                curr_pcc_add_onecase(i)=abs(corr(temp_add_onecase1',temp_add_onecase2'));
            end
            prob_add_onecase=curr_pcc_add_onecase/sum(curr_pcc_add_onecase);
            tmp_local_entropy_add_onecase=0;
            for i=1:e
                tmp_local_entropy_add_onecase=tmp_local_entropy_add_onecase-prob_add_onecase(i)*log(prob_add_onecase(i));
            end
            weighted_delta_entropy(na,l,s)=abs(tmp_local_entropy_add_onecase-tmp_local_entropy);
        end
        if mod(na,1000)==0
            na,l
        end
    end
    currtime=clock;
    l,etime(currtime,pretime)
    pretime=currtime;
end
%weighted_delta_entropy_e2=weighted_delta_entropy;
save weighted_delta_entropy.mat weighted_delta_entropy;
load('weighted_delta_entropy.mat');
close all;
aver_weighted_delta_entropy=zeros(psize(2),psize(3));
for l=1:psize(2)
    for s=1:patients_num(l)   %psize(3)
        tmp_delta_entropy=sort(weighted_delta_entropy(:,l,s),'descend');
        aver_weighted_delta_entropy(l,s)=mean(tmp_delta_entropy(1:floor(psize(1)*0.01)));
    end
    hold on;
    plot(l,aver_weighted_delta_entropy(l,1:patients_num(l)),'*');
    aver_aver(l)=mean(aver_weighted_delta_entropy(l,1:patients_num(l)));
    hold on;
end
figure(2);
plot([1:psize(2)],aver_aver);

fw=fopen('expression_to_deltaH.txt','w');
for i=1:psize(1)
    fprintf(fw,'%s',pipi{i});
    for s=1:patients_num(7)
        fprintf(fw,'\t%f',weighted_delta_entropy(i,7,s));
    end
    fprintf(fw,'\n');
end
fclose(fw);

fw=fopen('expression_to_deltaH_allStages.txt','w');
for i=1:psize(1)
    fprintf(fw,'%s',pipi{i});
    for l=1:psize(2)
        for s=1:patients_num(l)
            fprintf(fw,'\t%f',weighted_delta_entropy(i,l,s));
        end
    end
    fprintf(fw,'\n');
end
fclose(fw);


% [se,idx]=sort(weighted_delta_entropy(:,7,2),'descend');
% selected_deltaH_genes= pipi(idx(1:floor(psize(1)*0.05)));
% fid=fopen('selected_deltaH_genes.txt','w');
% %fprintf(fid,'%s\n',hline_patient_IDs{1}{1});
% for i=1:length(selected_deltaH_genes)
%     fprintf(fid,'%s\n',selected_deltaH_genes{i});
% end
% fclose(fid);
