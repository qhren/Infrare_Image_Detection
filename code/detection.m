function [re]=detection(f,count)
%����Ƶ�н�ȡ��һ֡ͼ�����ȶ���ͼ������
%ת��Ϊ�Ҷ�ͼ����룬��unit8����256*256�����ص�
% f=imread('./result/original_image/70.jpg');
[a,b]=size(f);
temp_Array=[ ];
%path_size��step_size��ѡȡ����ʵ�������趨
patch_size=80;
herizontal_size=6;
vertical_size=6;
PatchCount_w=ceil(((b-patch_size)/herizontal_size+1));
PatchCount_h=ceil((a-patch_size)/vertical_size+1);
element=zeros(patch_size*patch_size,PatchCount_w*PatchCount_h);
temp_Array=zeros(patch_size,patch_size);
for i=1:PatchCount_h  %��patch_size=80*80,herizontal_size=vertical_size=8ȡ�������patch_image
    for j=1:PatchCount_w
        if((j==PatchCount_w)&&(i==PatchCount_h))
         temp_Array(1:length((1+(i-1)*vertical_size):a),1:length((1+(j-1)*herizontal_size):b))=f((1+(i-1)*vertical_size):a,(1+(j-1)*herizontal_size):b);
        elseif((j==PatchCount_w)&&(i~=PatchCount_h))
         temp_Array=[f((1:patch_size)+(i-1)*vertical_size,(1+(j-1)*herizontal_size):b),zeros(patch_size,patch_size-length((1+(j-1)*herizontal_size):b))];
        elseif((i==PatchCount_h)&&(j~=PatchCount_w))
         temp_Array=[f((1+(i-1)*vertical_size):a,(1:patch_size)+(j-1)*herizontal_size);zeros(patch_size-length((1+(i-1)*vertical_size):a),patch_size)];
        else
         temp_Array=f((1:patch_size)+(i-1)*vertical_size,(1:patch_size)+(j-1)*herizontal_size);
        end
        element(:,(i-1)*(PatchCount_w)+j)=temp_Array(:);
    end
end
%%
clear temp_Array;
[M,N]=size(element);
%%
% solution via Accelerated Proximal Gradient
% ���ü����ݶ��㷨�������
% Initialization
B1=zeros(M,N);B2=zeros(M,N);T2=zeros(M,N);T1=zeros(M,N);A(2)=1;A(1)=1;%����ά����洢���Լ�С������
lamda=1/sqrt(max(a,b));
yita=0.9;%���ٵ�������
singular_values=svd(double(element),'econ');
singular_values_all=sort(singular_values);
u(2)=singular_values_all(length(singular_values)-1);
average_u=0.05*singular_values_all(length(singular_values)-3);
k=2;
%%
while abs(u(k)-average_u)>100  %ͨ�����Ƶ������ȿ��Լ�С���Ӷ�
    Yk1=B2+(A(k-1)-1)/A(k)*(B2-B1);
    Yk2=T2+(A(k-1)-1)/A(k)*(T2-T1);
    Gk=Yk1-1/2*(Yk1+Yk2-double(element));
    [U,S,V]=svd(Gk,'econ');%��ʡ������
    B3=U*wthresh(S,'s',u(k)/2)*V';
    temp1=B3;
    Gk2=Yk2-1/2*(Yk1+Yk2-double(element));
    T3=wthresh(Gk2,'s',lamda*u(k)/2);
    temp2=T3;
    A(k+1)=(1+sqrt(4*A(k)*A(k)+1))/2;
    u(k+1)=max(yita*u(k),average_u);
    k=k+1;
    B1=B2;
    B2=temp1;
    T1=T2;
    T2=temp2;
    disp(strcat('���ڽ��е�',num2str(k),'�ε�������'));
end
c=k-1;
B_1=B2;
T_1=T2;%�õ���B��T��patch_image��size��ͬ
%% ��������Ԫ����洢�ܽ��ʹ洢�����ǻ����޷����������ٶ�
%reconstruction of image
for i=1:PatchCount_h
    for j=1:PatchCount_w
    temp_image(:,:,(i-1)*PatchCount_w+j)=zeros(a,b);
    temp_image1(:,:,(i-1)*PatchCount_w+j)=zeros(a,b);
    tempr = reshape(T_1(:,(i-1)*PatchCount_w+j),[patch_size,patch_size]);
    tempr1=reshape(B_1(:,(i-1)*PatchCount_w+j),[patch_size,patch_size]);
    if((j==PatchCount_w)&&(i==PatchCount_h))
    temp_image((1+(i-1)*vertical_size):a,(1+(j-1)*herizontal_size):b,(i-1)*PatchCount_w+j)=tempr(1:length((1+(i-1)*vertical_size):a),1:length((1+(j-1)*herizontal_size):b));
    temp_image1((1+(i-1)*vertical_size):a,(1+(j-1)*herizontal_size):b,(i-1)*PatchCount_w+j)=tempr1(1:length((1+(i-1)*vertical_size):a),1:length((1+(j-1)*herizontal_size):b));
    elseif((j==PatchCount_w)&&(i~=PatchCount_h))
    temp_image((1:patch_size)+(i-1)*vertical_size,(1+(j-1)*herizontal_size):b,(i-1)*PatchCount_w+j)=tempr(1:patch_size,1:length((1+(j-1)*herizontal_size):b));
    temp_image1((1:patch_size)+(i-1)*vertical_size,(1+(j-1)*herizontal_size):b,(i-1)*PatchCount_w+j)=tempr1(1:patch_size,1:length((1+(j-1)*herizontal_size):b));
    elseif((j~=PatchCount_w)&&(i==PatchCount_h))
    temp_image((1+(i-1)*vertical_size):a,(1:patch_size)+(j-1)*herizontal_size,(i-1)*PatchCount_w+j)=tempr(1:length((1+(i-1)*vertical_size):a),1:patch_size);
    temp_image1((1+(i-1)*vertical_size):a,(1:patch_size)+(j-1)*herizontal_size,(i-1)*PatchCount_w+j)=tempr1(1:length((1+(i-1)*vertical_size):a),1:patch_size);
    else
    temp_image((1:patch_size)+(i-1)*vertical_size,(1:patch_size)+(j-1)*herizontal_size,(i-1)*PatchCount_w+j)=tempr(1:patch_size,1:patch_size);
    temp_image1((1:patch_size)+(i-1)*vertical_size,(1:patch_size)+(j-1)*herizontal_size,(i-1)*PatchCount_w+j)=tempr1(1:patch_size,1:patch_size);
    end
    end
end
%% ��ֵ�˲�
re_B=zeros(a,b);
re_T=zeros(a,b);
result_T=zeros(a,b,N);
result_B=zeros(a,b,N);
for i=1:a
    for j=1:b
       for  k=1:N
          result_T(i,j,k)=temp_image(i,j,k);
          result_B(i,j,k)=temp_image1(i,j,k);
       end
       if(isempty(find(result_T(i,j,:)~=0, 1)))
       temp2=0;
       re_T(i,j)=0;
       else
       temp2=result_T(i,j,:);
       position2=find(temp2~=0);
       re_T(i,j)=median(reshape(temp2(temp2~=0),[1,length(position2)]));
       end
       if(isempty(find(result_B(i,j,:)~=0, 1)))
       temp1=0;
       re_B(i,j)=0;
       else
       temp1=result_B(i,j,:);
       position1=find(temp1~=0);
       re_B(i,j)=median(reshape(temp1(temp1~=0),[1,length(position1)]));
      end
    end
end
%% ���޴���
k=1;
[row,col]=find(re_T);
for i=1:length(row)
    if (re_T(row(i),col(i))>max(mean(re_T(re_T~=0)),mean2(re_T)+c*std2(re_T)))% attention to bright targets
        re_T(row(i),col(i))=255;
         row_new(k)=row(i);
         col_new(k)=col(i);
         k=k+1;
    else
        re_T(row(i),col(i))=0;
    end
end
if(k==1)
    row_new=[];
    col_new=[];
end
%%
% �����о������о��󣬽�Ŀ���λ����ԭʼ��ͼ���б��
%(�����������Ŀ�����ܳ�����ͼ��߽磬��Ƕ�����������������ã�
% if(~isempty(row_new))
% for i=1:length(row_new)
%     f(row_new(i),(col_new(i)-7):(col_new(i)-5))=256*ones(1,3);
%     f(row_new(i),(col_new(i)+5):(col_new(i)+7))=256*ones(1,3);
%     f((row_new(i)-7):(row_new(i)-5),col_new(i))=256*ones(1,3);
%     f((row_new(i)+5):(row_new(i)+7),col_new(i))=256*ones(1,3);    
% end
% end
% imshow(f);
re=re_T;
% imwrite(f,strcat('D:\Infraed_image\avi_detection\result\original_mark\',num2str(count),'.jpg'),'jpg');%��������� д���ļ�
end



        
