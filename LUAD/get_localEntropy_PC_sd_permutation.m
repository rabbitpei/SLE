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
D=randperm(size(tempcontrol,2));
% weighted_delta_entropy=zeros(total_node_num,psize(2),psize(3));
weighted_delta_sd=zeros(total_node_num,psize(2),psize(3));
for l=1:psize(2)
    for na=1:total_node_num
        edges_list=[];
        neighbour_num=length(adjacent_network{na})-1;
        center=str2num(adjacent_network{na}{1});
        sd_control=std(tempcontrol(center,D));
        for s=1:patients_num(l)
            tmp_sd_add_onecase=std([tempcontrol(center,D),...
                reshape(tempcase(center,l,s),1,1)]);
            weighted_delta_sd(na,l,s)=abs(tmp_sd_add_onecase-sd_control);
            weighted_edges(na,l,s)=log(neighbour_num);
        end
    end
    
    currtime=clock;
    l,etime(currtime,pretime)
    pretime=currtime;
end
save weighted_delta_sd;
save weighted_edges;
load weighted_delta_entropy;

for i=1:total_node_num
    for l=1:psize(2)
        for s=1:patients_num(l)
            complex_index(i,l,s)=weighted_delta_entropy(i,l,s);  %*weighted_delta_sd(i,l,s);
            complex_index2(i,l,s)=weighted_delta_sd(i,l,s)*log(weighted_edges(i,l,s));
            complex_index3(i,l,s)=weighted_delta_entropy(i,l,s)*weighted_delta_sd(i,l,s)/log(weighted_edges(i,l,s));  %/weighted_edges(i,l,s);
        end
    end
end

save complex_index3;
aver_comidx=zeros(psize(2),psize(3));
for l=1:psize(2)
    for s=1:patients_num(l)
        tmp_comidx=sort(complex_index(:,l,s),'descend');
        aver_comidx(l,s)=mean(tmp_comidx(1:floor(psize(1)*0.01)));
        
        tmp_comidx2=sort(complex_index2(:,l,s),'descend');
        aver_comidx2(l,s)=mean(tmp_comidx2(1:floor(psize(1)*0.01)));
        
        tmp_comidx3=squeeze(complex_index3(:,l,s));
        tmp_comidx3(find(isnan(tmp_comidx3)))=0;
        %tmp_comidx3=sort(complex_index3(:,l,s),'descend');
        aver_comidx3(l,s)=mean(tmp_comidx3(1:floor(psize(1)*0.01)));
    end
    aver_aver(l)=mean(aver_comidx(l,1:patients_num(l)));
    aver_aver2(l)=mean(aver_comidx2(l,1:patients_num(l)));
    aver_aver3(l)=mean(aver_comidx3(l,1:patients_num(l)));
    aver_aver3_2(l)=aver_aver(l)*aver_aver2(l);
%     aver_sd1(l)=std(aver_comidx(l,1:patients_num(l)));
%     aver_sd2(l)=std(aver_comidx2(l,1:patients_num(l)));
%     aver_sd3(l)=std(aver_comidx3(l,1:patients_num(l)));
    %aver_aver3(l)=aver_aver(l)*aver_aver2(l);
end

close all;
plot([1:psize(2)],aver_aver);
figure(2)
plot([1:psize(2)],aver_aver2);
figure(3)
plot([1:psize(2)],aver_aver3);
figure(4)
plot([1:psize(2)],aver_aver3_2);
for l=7:7
    tmp_comidx=mean(reshape(complex_index3(:,l,1:patients_num(l)),psize(1),patients_num(l)),2);
    tmp_comidx(find(isnan(tmp_comidx)))=0;
    [se,idx]=sort(tmp_comidx,'descend');
    %aver_ci(l)=mean(tmp_comidx(idx(1:floor(psize(1)*0.01)))); 
    top_comidx_genes= pipi(idx(1:floor(psize(1)*0.03)));    
end

fid=fopen('selected_comidx_genes_0.01.txt','w');
for i=1:length(top_comidx_genes)
    fprintf(fid,'%s\n',top_comidx_genes{i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fid=fopen('selected_comidx_genes_data_0.03.txt','w');
fid2=fopen('selected_comidx_genes_data_original_0.03.txt','w');
for l=1:8
    tmp_comidx_each(:,l)=mean(reshape(complex_index3(:,l,1:patients_num(l)),psize(1),patients_num(l)),2);
    tmp_comidx_each(find(isnan(tmp_comidx)),l)=0;
    %tmp_comidx_each(:,l)=mean(reshape(complex_index3(:,2,1:patients_num(l)),psize(1),patients_num(l)),2);
end
tmp_a=zeros(length(top_comidx_genes),8);
for i=1:length(top_comidx_genes)
    mina=min(tmp_comidx_each(idx(i),:));
    maxa=max(tmp_comidx_each(idx(i),:));
    tmp_a(i,:)=(tmp_comidx_each(idx(i),:)-mina)/(maxa-mina);
    tmp=zscore(tmp_comidx_each(idx(i),:));
    fprintf(fid2,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',top_comidx_genes{i},tmp(2),tmp(3),tmp(4),tmp(5),tmp(6),tmp(7),tmp(8));
end
sl=0;
for l=2:2
    for i=1:length(top_comidx_genes)
        if tmp_a(i,l)<tmp_a(i,7) && tmp_a(i,7)>0.5
            sl=sl+1;
            fprintf(fid,'%s\t%f\t%f\n',top_comidx_genes{i},tmp_a(i,2),tmp_a(i,7));
        end
    end
%     if sl>100
%         sl,l
%         break;
%     end
end
fclose(fid2);
fclose(fid);


gene_idx_count=zeros(psize(1),1);
for l=7:7
    for j=1:patients_num(l)
        fid=fopen(['selected_genes_0.03_sample',num2str(j),'.txt'],'w');
        tmp=squeeze(complex_index3(:,l,j));
        tmp(find(isnan(tmp)))=zeros(1,length(tmp(find(isnan(tmp)))));
        [se,idx]=sort(tmp,'descend');
        top_genes=pipi(idx(1:floor(psize(1)*0.03)));
        for i=1:floor(psize(1)*0.03)
            fprintf(fid,'%s\n',top_genes{i});
        end
        tmp_idx=idx(1:floor(psize(1)*0.03));
        gene_idx_count(tmp_idx)=gene_idx_count(tmp_idx)+1;
        
        
        fclose(fid);
    end
    
end