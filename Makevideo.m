function Makevideo(Path,Outpath)
%将指定路径Path下的图片制作成视频,命名为videoname，E:\Luoxy\clump_core\result4
% 保存在OutPath路径下，参数均为字符串格式G:\远端增亮2013\201301010043_0131_B
% Path='G:\hsos\2014\01jpg';
% 获取视频的名称
videoname=split(Path,'\');
videoname=[videoname{end-1} '_' videoname{end}];

fdt1=dir(Path);
% 保证非空
if length(fdt1)>3
    FileStyle=split(fdt1(3).name,'.');
    FileStyle=['*.' FileStyle{end}];    % 获取文件的类型
   
    if strcmpi(FileStyle,'*.png')
        aviobj = VideoWriter([Outpath '\' videoname '.avi']);
        aviobj.FrameRate =10; %代表1秒播放2帧图片
        open(aviobj)
        for j=3:size(fdt1,1)
            fitsname=fdt1(j).name;
            frame = imread([Path '\' fitsname]);
            frame =mat2gray(frame);
%             frame=GuiYiHua(frame);
%             frame=frame(1:1900,1:1900);
            writeVideo(aviobj,frame);
        end
        close(aviobj)
    end
    
    if strcmpi(FileStyle,'*.fits')
        aviobj = VideoWriter([Outpath '\' videoname '.avi']);
        aviobj.FrameRate =10; %代表1秒播放2帧图片
        open(aviobj)
        for j=3:size(fdt1,1)
            fitsname=fdt1(j).name;
            frame = fitsread([Path '\' fitsname]);
            frame =mat2gray(frame(:,:,1));
            frame=GuiYiHua(frame);
            writeVideo(aviobj,frame);
        end
        close(aviobj)
    end
end



% numzeros= 2;    %视频名字的长度
% nz = strcat('%0',num2str(numzeros),'d');
%  id=sprintf(nz,k);
% % mkdir('data1')
% for k = 1 : numFrames% 读取前15帧
%     frame = read(Moveobj,k);%读取第几帧
%     id=sprintf(nz,k);
%     imwrite(frame,strcat('data1\',id,'.jpg'),'jpg');% 保存帧
% end


