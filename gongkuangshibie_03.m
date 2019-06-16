function []=errerere();
close all;
global data sample scale;
if numel(data)==0;
    %读取xls数据
    file1='标准工况数据.xls';
    file2='实车行驶工况数据.xlsx';
    range1='A3:AN15833';
    range2='A2:B497';
    [data,TXT,RAW]=xlsread(file1,1,range1);
    [sample,TXT,RAW]=xlsread(file2,1,range2);
end

n=size(data,2)/2;
%提取特征
feature=[];
figure(1);
for i=1:n;
    tt=data(:,i*2-1);
    vv=data(:,i*2);
    [featurei,tt,vv]=get_feature(tt,vv);
    subplot(5,4,i);
    plot(tt,vv,'b.-');title(strcat('mode-',num2str(i)))
    if isnan(vv(end));
        error;
    end
    feature=[feature;featurei(:)'];
end
disp('原始特征参数(17类)')
disp(feature),

%对不同属性赋予不同的权重
w=sqrt([0.44,0.10,0.34,0.04,0.08]);%权重，参考(b-7)
%w=[2,2,2,1,1];
% w=[2,0.01,0.01,0.01,0.01];
global ceshi_index;
ceshi_index=1;
w=0.001*ones(1,5);w(ceshi_index)=1;

%归一化
scale=(max(feature)-min(feature));
scale=max(1e-3,scale./w);%公式(b-7)
feature_scale=feature./(ones(n,1)*scale);

%聚类分析
%dis=pdist(feature_scale(:,1:3),'euclidean);
dis=pdist(feature_scale,'euclidean');
link=linkage(dis);
figure;dendrogram(link);xlabel('循环工况序号');ylabel('综合聚类距离尺度')
c=cluster(link,'maxclust',5);
li=feature;%li=feature_scale;
figure;scatter3(li(:,1),li(:,2),li(:,3),200,c,'filled');%hold on;for i=1:numel(c);annotation('textbox',[li(i,:)],);end
name={'max(V)','mean(V)','std(V)'};title('聚类结果');xlabel(name{1});ylabel(name{2});zlabel(name{3});
if 1
    %从cluster聚类结果中得出聚类
    mode={};
    uc=unique(c);
    for i=1:numel(uc);
        mode{i}=find(c==uc(i));
    end
end

%提取工况特征值矩阵 公式(b-2)
Y=[];Y1=[];
for i=1:numel(mode);
    index=mode{i};
    feature_scale_i=feature_scale(index,:);
    feature_i=feature(index,:);
    if numel(index)>1;yi=mean(feature_scale_i);
        yi1=mean(feature_i);
    else
        yi=feature_scale_i;
        yi1=feature_i;
    end
    Y=[Y,yi(:)];
    Y1=[Y1,yi1(:)];
end

%根据测试属性(ceshi_index)大小给每一类赋予标签
[~,li]=sort(Y1(ceshi_index,:));
label=[];for i=1:numel(li);label(i)=find(li==i);end;

%显示聚类结果
disp('聚类后的特征参数(5类)');
disp('类别标签');
disp(label);
disp('特征参数值');
disp(Y1)
disp('尺度缩放后的特征参数(5类)');
disp(Y)
disp('分类情况');for i=1:numel(mode);disp(strcat('第',num2str(i),'类包含'));disp(mode{i}');end;

global Yoriginal;Yoriginal=Y1;
%计算Y 公式(b-5)，由于之前已经归一化了，这里不再需要归一化，直接用就行
%从实际工况数据中提取特征
[~,tt,vv]=get_feature(sample(:,1),sample(:,2));%从实际工况数据中提取特征%
%  k=3;[~,tt,vv]=get_feature(data(:,k*2-1),data(:,k*2));%从标准工况数据中提取特征
%vv=pinghua(vv);
%计算地鼠度S 公式(b-6)
cm=[5,10,15,20];
figure;
global ceshi;ceshi=[];
p=0;nfeature=0;
subplot(numel(cm)+(nfeature)-p+1,1,1);plot(tt,vv,'linewidth',3);ylabel('车速/km/h');set(gca,'xtick',[]);
name={'最大速度特征参数','平均速度特征','速度标准差','平均正加速度特征','平均负加速度特征'};
for i=1:numel(cm);
    m_horizon=cm(i);%时间窗长度
    [cmode,t1,xx,aa]=mode_identify(tt,vv,m_horizon,scale,Y,label);%识别模式
    subplot(numel(cm)+(nfeature)-p+1,1,i+1);plot(t1,cmode,'o-')
    grid on
    axis([0,Inf,1,numel(mode)]);
    set(gca,'Ytick',[1:numel(mode)]);
    set(gca,'xtick',[]);
    %if i==numel(cm);xlabel('时间/s');end
    ylabel(strcat('识别(',num2str(cm(i)*10),'s)'));%title(strcat('时间窗长度=',num2str(cm(i)),'s'));
end
xlabel('时间/s');
%figure;plot(tt,aa)

1;
%figure;plot(ceshi,'linewidth',5);legend('当时最大速度','匹配工况最大速度特征参数');
%figure;plot(t1,xx(:,ceshi_index));

function [cmode,t1,xx,aa]=mode_identify(tt,vv,m_horizon,scale,Y,label)
%识别模式
global Yoriginal;
global ceshi;
global ceshi_index;
n=numel(tt);
uu=[];t1=[];cmode=[];
xx=[];
aa=jiasudu(tt,vv);
for i=1:n-m_horizon;
    index=i:i+m_horizon;
    ti=tt(index);
    vi=vv(index);
    ai=aa(index);
    xi=get_feature(ti,vi,ai);%提取特征参数
    if sum(isnan(xi))>0;
        1;
    end
    ui=calU(xi,Y,scale);%计算地鼠度
    uu=[uu;ui(:)'];%地鼠度
    [~,k]=max(ui);
    if 1
        %这是测试用的代码，可以忽略
        ceshi=[ceshi;[xi(ceshi_index),Yoriginal(ceshi_index,k)]];
        if (xi(ceshi_index)<57)&&(Yoriginal(ceshi_index,k)>70);
            1;
        end
    end
    cmode=[cmode;label(k)];%label(k)是大类编号
    t1=[t1;ti(end)];%时间
    xx=[xx;xi(:)'];
end
ceshi;

function U=calU(x,Y,scale)
%计算x的地鼠度
%Y是标准工况的属性矩阵
%scale是缩放尺度
%采用欧氏距离
x=x(:)'./scale;%缩放实际工况参数
S=squareform(pdist([x',Y]', 'euclidean'));%计算实际工况和标准工况之间的距离
ep=1e-8;U=1./(ep+S(1,2:end));%地鼠度是距离的反函数
%根据距离计算相对地鼠度
U=U/sum(U);%归一化

function out=distance(x,y,w,p)
%公式（9）中的距离定义
out=(w.*(x-y).^p).^(1/p);
function [out,tt,vv]=get_feature(tt,vv,a)
k=find(isnan(tt));
if numel(k)>0;
    k=k(1)-1;
else
    k=numel(tt);
end
tt=tt(1:k,:);
vv=vv(1:k,:);

vmax=max(vv);
vav=mean(vv);
vstd=std(vv);
%加速度提取
ep=0.0;
if 1
    if nargin<3;
        a=jiasudu(tt,vv);
    end
    ap=mean(a(a>=ep));an=-mean(a(a<=-ep));
end
if isnan(ap);
    ap=0;
end
if isnan(an);
    an=0;
end
%apstd=std(a(a>ep));anstd=std(a(a<-ep));
out=[vmax,vav,vstd,ap,an];
function vv1=pinghua(vv,c,k);
n=numel(vv);
vv1=vv;
if nargin<2;
c=5;
end
if nargin<3;
    k=0.8;
end
for i=1:n-2*c-1;
   ind=i:i+2*c+1;
    vv1(i+c)=(sum(vv(ind))-vv(i+c))/2/c*k+(1-k)*vv(i+c);
end
function aa=jiasudu(tt,vv);
vv1=pinghua(vv,1,0.5);
n=numel(vv);
c=3;
aa=vv*0;
for i=1:n-2*c-1;
    aa(i+c)=(vv(i+c*2+1)-vv(i))/(tt(i+c*2+1)-tt(i));
end
aa=pinghua(aa,1,0.5);