
clear;
clc;
close all;
%load alldata;
load T1;
%[T1] = xlsread('data32_filled.xlsx');
load genename;
profile=T1;
% profile=zscore(profile');
% profile=profile';

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
        center=str2num(adjacent_network{na}{1});
        neighbour_num=length(adjacent_network{na});
        for s=1:psize(3)
            sd_control=std(reshape(tempcontrol(center,D,s),5,1));
            tmp_sd_add_onecase=std([reshape(tempcontrol(center,D,s),1,5),...
                reshape(tempcase(center,l,s),1,1)]);
            weighted_delta_sd(na,l,s)=abs(tmp_sd_add_onecase-sd_control)*log(tmp_sd_add_onecase/sd_control);
        end 
        
    end
    
    currtime=clock;
    l,etime(currtime,pretime)
    pretime=currtime;
end

weighted_delta_sd2=weighted_delta_sd;
save weighted_delta_sd weighted_delta_sd;
load weighted_delta_entropy;
ctrl_patients_idx=[2,3,4,9,11,14,16,17];
case_patients_idx=[1,5,6,7,8,10,12,13,15];
close all;
figure(1);
for s=2:length(case_patients_idx)
    ss=case_patients_idx(s);
    for l=1:psize(2)
        tmp_delta_sd=sort(weighted_delta_sd(:,l,ss),'descend');
        aver_weighted_delta_sd(l)=mean(tmp_delta_sd(1:floor(total_node_num*0.01)));
    end
    plot([1:psize(2)],aver_weighted_delta_sd,'r','LineWidth',3);
    hold on;
end
%title('Symptomatic samples SD');
%figure(2)
for s=1:length(ctrl_patients_idx)
    for l=1:psize(2)
        ss=ctrl_patients_idx(s);
        tmp_delta_sd=sort(weighted_delta_sd(:,l,ss),'descend');
        aver_weighted_delta_sd(l)=mean(tmp_delta_sd(1:floor(total_node_num*0.01)));
    end
    plot([1:psize(2)],aver_weighted_delta_sd,'b','LineWidth',3);
    hold on;
end

for i=1:total_node_num
    for l=1:psize(2)
        for s=1:psize(3)
            complex_index(i,l,s)=weighted_delta_sd(i,l,s)*weighted_delta_entropy(i,l,s);
        end
    end
end
save complex_index;

figure('NumberTitle', 'off', 'Name', 'Global Complex-index curves');
CI=zeros(psize(2),psize(3));
for s=1:length(case_patients_idx)
    ss=case_patients_idx(s);
    for l=1:psize(2)
        tmp_com_idx=sort(complex_index(:,l,ss),'descend');
        aver_com_idx(l)=mean(tmp_com_idx(1:floor(total_node_num*0.01)));
    end
    CI(:,ss)=aver_com_idx(1:psize(2));
end

for s=1:length(ctrl_patients_idx)
    ss=ctrl_patients_idx(s);
    for l=1:psize(2)
        tmp_com_idx=sort(complex_index(:,l,ss),'descend');
        aver_com_idx(l)=mean(tmp_com_idx(1:floor(total_node_num*0.01)));
    end
    CI(:,ss)=aver_com_idx(1:psize(2));
end
xlim([1,psize(2)]);
plot(CI(:,case_patients_idx(1)),'r','LineWidth',3);
hold on;
plot(CI(:,ctrl_patients_idx(1)),'b','LineWidth',3);
legend('Sx','Asx');
plot(CI(:,case_patients_idx),'r','LineWidth',3);
plot(CI(:,ctrl_patients_idx),'b','LineWidth',3);

%%%%%%%%%%%%%%%%%%%%%%%  landscape individual %%%%%%%%%%%%%%%
figure('NumberTitle', 'off', 'Name', 'Individual landscape: Symptom');
tipping_point=[6,8,7,8,9,11,8,10,12];
for s=1:9
    subplot(3,3,s);
    [tmp_com_idx,idx]=sort(complex_index(:,tipping_point(s),s),'descend');
    a1=squeeze(complex_index(idx(1:2000),:,s));
    
    Y=pdist(a1(:,s),@cofun2);
    Z=linkage(Y);
    [H,ter,perm] = dendrogram(Z, 0, 'orientation', 'left', ...
        'colorthreshold', 'default');
    a2=a1(perm,:);
    surf([1:16],[1:200],a1(1:200,:));
    shading interp;
end
%%%%%%%%%%%%%%%%%%%%% landscape ############
figure('NumberTitle', 'off', 'Name', 'Global CI landscape');
LCI(:,1:length(case_patients_idx))=CI(:,case_patients_idx);
LCI(:,length(case_patients_idx)+1:length(case_patients_idx)+length(ctrl_patients_idx))=CI(:,ctrl_patients_idx);
surf([1:size(CI,2)],[1:size(CI,1)],LCI);
shading interp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tc=[8,8,8,8,11,12,14,15,15];
sx=CI(:,case_patients_idx);
asx=CI(:,ctrl_patients_idx);
va=ceil(max(max(CI)));
n1=0;
correct1=[];
stl=ceil(va/4000);
for l=0:stl:stl*4000
    for i=1:length(case_patients_idx)
        for j=1:psize(2)
            if sx(j,i)>l&&j<tc(i)
                n1=n1+1;
                break
            end
        end
    end
    correct1=[correct1,n1];
    n1=0;
end
correct1(2,:)=[0:stl:stl*4000];

n2=0;
correct2=[];
for l=0:stl:stl*4000
    for i=1:8
        for j=1:psize(2)
            if asx(j,i)>l%I2??timepoint???
                n2=n2+1;
                break;
            end
        end
    end
    correct2=[correct2,n2];%???
    n2=0;
end

y=correct1(1,:)/9;%???
x=correct2/8;%???
% x(61)=0;
% y(61)=0;
figure1=figure('color',[1 1 1],'NumberTitle', 'off', 'Name', 'ROC');
ylim([0,1])
plot(x,y,'LineWidth',3)
% legend('ROC')
s=0;
for i=1:4000
    s1=1/2*abs(y(i+1)+y(i))*abs(x(i+1)-x(i));
    s=s+s1;
end


figure('NumberTitle', 'off', 'Name', 'Individual Complex-index curves: Symptom');
%%%???
ac=correct1(1,:)-correct2(1,:);
v1=find(ac==max(ac));
v=stl*v1(1);
%v=300;
for j=1:9
    subplot(3,3,j)
    plot(1:tn,sx(:,j),'k-','LineWidth',2);
    hold on
    plot(tc(j),sx(tc(j),j),'mo','markersize',12,'LineWidth',2);
    
    V = find(sx(:,j)>=v+1);%%%??????????????????
    if size(V,1)>0
        plot(V(1),sx(V(1),j),'bP','markersize',tn,'LineWidth',2);
    end
     
    xlim([1,tn]);
    xlabel('Time');
    ylabel('Early-warning Signal');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('NumberTitle', 'off', 'Name', 'Individual Complex-index curves: Asymptom');
%v=300;
for j=1:length(ctrl_patients_idx)
    subplot(3,3,j)
    plot(1:tn,asx(:,j),'k-','LineWidth',2);
    hold on
    plot(tc(j),asx(tc(j),j),'mo','markersize',12,'LineWidth',2);
    
    V = find(asx(:,j)>=v+1);%%%??????????????????
    if size(V,1)>0
        plot(V(1),asx(V(1),j),'bP','markersize',tn,'LineWidth',2);
    end
     
    xlim([1,tn]);
    xlabel('Time');
    ylabel('Early-warning Signal');
end





%%%%%%%%% selected marker genes %%%%%%%%%%%%%%%%%%%%%%%%%%%%
selected_deltaH_genes=[];
fid=fopen('selected_flu_genes_0.05.txt','w');
selected_deltaH_genes=containers.Map;
predicted.point=[6,8,7,8,9,11,8,10,12];
for s=1:length(case_patients_idx)
    ss=case_patients_idx(s);
    [tmp_com_idx,idx]=sort(complex_index(:,predicted.point(s),ss),'descend');
    %aver_com_idx(l)=mean(tmp_com_idx(1:floor(total_node_num*0.01)));
    
    tmp_genes=symbols(idx(1:floor(total_node_num*0.05)));
    for i=1:length(tmp_genes)
        if isKey(selected_deltaH_genes,tmp_genes(i))==0
            selected_deltaH_genes(tmp_genes{i})=1;
        else
            selected_deltaH_genes(tmp_genes{i})=selected_deltaH_genes(tmp_genes{i})+1;
        end
    end
end
genes=selected_deltaH_genes.keys();
for i=1:length(genes)
    fprintf(fid,'%s\t%d\n',genes{i},selected_deltaH_genes(genes{i}));
end
fclose(fid);
