Obj = VideoReader('yuanshi.avi');%��ȡ��Ƶ  
vidWidth = Obj.Width;  
vidHeight = Obj.Height;  
video = struct('frame',zeros(vidHeight,vidWidth,3,'uint8'));  
k = 1;  
while hasFrame(Obj)  
    video(k).frame = rgb2gray(readFrame(Obj));  
   imwrite(video(k).frame,strcat('D:\Infraed_image\avi_detection\result\original_image\',num2str(k),'.jpg'),'jpg');% ����ÿһ֡
    k = k+1;
end
for count=51:length(video)
    re=detection(video(count).frame,count);
    imwrite(mat2gray(re),strcat('D:\Infraed_image\avi_detection\result\Seperate_target(80-80-6)\',num2str(count),'.jpg'),'jpg');% ����ÿһ֡
    count
end
detectobj = VideoWriter('process_out.avi');%Ŀ����Ƶ·��������
writerObj.FrameRate = 30;
open(detectobj);
for i=1:86%��֡������Ƶ��
    fname=strcat('D:\Infraed_image\avi_detection\result\Seperate_target(80-80-6)\',num2str(i),'.jpg');
    frame = imread(fname);
    writeVideo(detectobj,frame);
end 
close(detectobj);%�ر���Ƶ����


