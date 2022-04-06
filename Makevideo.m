function Makevideo(Path,Outpath)
%��ָ��·��Path�µ�ͼƬ��������Ƶ,����Ϊvideoname��E:\Luoxy\clump_core\result4
% ������OutPath·���£�������Ϊ�ַ�����ʽG:\Զ������2013\201301010043_0131_B
% Path='G:\hsos\2014\01jpg';
% ��ȡ��Ƶ������
videoname=split(Path,'\');
videoname=[videoname{end-1} '_' videoname{end}];

fdt1=dir(Path);
% ��֤�ǿ�
if length(fdt1)>3
    FileStyle=split(fdt1(3).name,'.');
    FileStyle=['*.' FileStyle{end}];    % ��ȡ�ļ�������
   
    if strcmpi(FileStyle,'*.png')
        aviobj = VideoWriter([Outpath '\' videoname '.avi']);
        aviobj.FrameRate =10; %����1�벥��2֡ͼƬ
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
        aviobj.FrameRate =10; %����1�벥��2֡ͼƬ
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



% numzeros= 2;    %��Ƶ���ֵĳ���
% nz = strcat('%0',num2str(numzeros),'d');
%  id=sprintf(nz,k);
% % mkdir('data1')
% for k = 1 : numFrames% ��ȡǰ15֡
%     frame = read(Moveobj,k);%��ȡ�ڼ�֡
%     id=sprintf(nz,k);
%     imwrite(frame,strcat('data1\',id,'.jpg'),'jpg');% ����֡
% end


