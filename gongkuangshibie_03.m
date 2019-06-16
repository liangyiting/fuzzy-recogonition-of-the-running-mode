function []=errerere();
close all;
global data sample scale;
if numel(data)==0;
    %��ȡxls����
    file1='��׼��������.xls';
    file2='ʵ����ʻ��������.xlsx';
    range1='A3:AN15833';
    range2='A2:B497';
    [data,TXT,RAW]=xlsread(file1,1,range1);
    [sample,TXT,RAW]=xlsread(file2,1,range2);
end

n=size(data,2)/2;
%��ȡ����
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
disp('ԭʼ��������(17��)')
disp(feature),

%�Բ�ͬ���Ը��費ͬ��Ȩ��
w=sqrt([0.44,0.10,0.34,0.04,0.08]);%Ȩ�أ��ο�(b-7)
%w=[2,2,2,1,1];
% w=[2,0.01,0.01,0.01,0.01];
global ceshi_index;
ceshi_index=1;
w=0.001*ones(1,5);w(ceshi_index)=1;

%��һ��
scale=(max(feature)-min(feature));
scale=max(1e-3,scale./w);%��ʽ(b-7)
feature_scale=feature./(ones(n,1)*scale);

%�������
%dis=pdist(feature_scale(:,1:3),'euclidean);
dis=pdist(feature_scale,'euclidean');
link=linkage(dis);
figure;dendrogram(link);xlabel('ѭ���������');ylabel('�ۺϾ������߶�')
c=cluster(link,'maxclust',5);
li=feature;%li=feature_scale;
figure;scatter3(li(:,1),li(:,2),li(:,3),200,c,'filled');%hold on;for i=1:numel(c);annotation('textbox',[li(i,:)],);end
name={'max(V)','mean(V)','std(V)'};title('������');xlabel(name{1});ylabel(name{2});zlabel(name{3});
if 1
    %��cluster�������еó�����
    mode={};
    uc=unique(c);
    for i=1:numel(uc);
        mode{i}=find(c==uc(i));
    end
end

%��ȡ��������ֵ���� ��ʽ(b-2)
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

%���ݲ�������(ceshi_index)��С��ÿһ�ำ���ǩ
[~,li]=sort(Y1(ceshi_index,:));
label=[];for i=1:numel(li);label(i)=find(li==i);end;

%��ʾ������
disp('��������������(5��)');
disp('����ǩ');
disp(label);
disp('��������ֵ');
disp(Y1)
disp('�߶����ź����������(5��)');
disp(Y)
disp('�������');for i=1:numel(mode);disp(strcat('��',num2str(i),'�����'));disp(mode{i}');end;

global Yoriginal;Yoriginal=Y1;
%����Y ��ʽ(b-5)������֮ǰ�Ѿ���һ���ˣ����ﲻ����Ҫ��һ����ֱ���þ���
%��ʵ�ʹ�����������ȡ����
[~,tt,vv]=get_feature(sample(:,1),sample(:,2));%��ʵ�ʹ�����������ȡ����%
%  k=3;[~,tt,vv]=get_feature(data(:,k*2-1),data(:,k*2));%�ӱ�׼������������ȡ����
%vv=pinghua(vv);
%��������S ��ʽ(b-6)
cm=[5,10,15,20];
figure;
global ceshi;ceshi=[];
p=0;nfeature=0;
subplot(numel(cm)+(nfeature)-p+1,1,1);plot(tt,vv,'linewidth',3);ylabel('����/km/h');set(gca,'xtick',[]);
name={'����ٶ���������','ƽ���ٶ�����','�ٶȱ�׼��','ƽ�������ٶ�����','ƽ�������ٶ�����'};
for i=1:numel(cm);
    m_horizon=cm(i);%ʱ�䴰����
    [cmode,t1,xx,aa]=mode_identify(tt,vv,m_horizon,scale,Y,label);%ʶ��ģʽ
    subplot(numel(cm)+(nfeature)-p+1,1,i+1);plot(t1,cmode,'o-')
    grid on
    axis([0,Inf,1,numel(mode)]);
    set(gca,'Ytick',[1:numel(mode)]);
    set(gca,'xtick',[]);
    %if i==numel(cm);xlabel('ʱ��/s');end
    ylabel(strcat('ʶ��(',num2str(cm(i)*10),'s)'));%title(strcat('ʱ�䴰����=',num2str(cm(i)),'s'));
end
xlabel('ʱ��/s');
%figure;plot(tt,aa)

1;
%figure;plot(ceshi,'linewidth',5);legend('��ʱ����ٶ�','ƥ�乤������ٶ���������');
%figure;plot(t1,xx(:,ceshi_index));

function [cmode,t1,xx,aa]=mode_identify(tt,vv,m_horizon,scale,Y,label)
%ʶ��ģʽ
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
    xi=get_feature(ti,vi,ai);%��ȡ��������
    if sum(isnan(xi))>0;
        1;
    end
    ui=calU(xi,Y,scale);%��������
    uu=[uu;ui(:)'];%�����
    [~,k]=max(ui);
    if 1
        %���ǲ����õĴ��룬���Ժ���
        ceshi=[ceshi;[xi(ceshi_index),Yoriginal(ceshi_index,k)]];
        if (xi(ceshi_index)<57)&&(Yoriginal(ceshi_index,k)>70);
            1;
        end
    end
    cmode=[cmode;label(k)];%label(k)�Ǵ�����
    t1=[t1;ti(end)];%ʱ��
    xx=[xx;xi(:)'];
end
ceshi;

function U=calU(x,Y,scale)
%����x�ĵ����
%Y�Ǳ�׼���������Ծ���
%scale�����ų߶�
%����ŷ�Ͼ���
x=x(:)'./scale;%����ʵ�ʹ�������
S=squareform(pdist([x',Y]', 'euclidean'));%����ʵ�ʹ����ͱ�׼����֮��ľ���
ep=1e-8;U=1./(ep+S(1,2:end));%������Ǿ���ķ�����
%���ݾ��������Ե����
U=U/sum(U);%��һ��

function out=distance(x,y,w,p)
%��ʽ��9���еľ��붨��
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
%���ٶ���ȡ
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