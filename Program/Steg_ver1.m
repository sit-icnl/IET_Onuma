
%��ҁF����C�m�@'20�� af16023
%               '22�� ma20017

%�X�y�N�g���g�U�@��(k,n)臒l�閧���U�@��p�����X�e�K�m�O���t�B
%�Í����T�C�g�Fhttp://point-at-infinity.org/ssss/demo.html




clc;
clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�@��v�p�����[�^�ݒ�@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ߍ��݊J�n�ʒu
ss = 5;  %MAX64-Leg+1

%DCT�W���̒��o��(���ɉ����Ė��ߍ��ޏ��ʂ��ő���܂œ���Ȃ��ƃG���[)
Leg = 2; %64
%Leg = 4; %128
%Leg = 8; %256
%Leg = 16; %512
%Leg = 32; %1024
%Leg = 50; %1600

%M'�n��̃Q�C��
%G = 4.159;   %�Q�C��
%G = 7.699;   %�Q�C��
%G = 11.036;  %�Q�C��
%G = 19.127;  %�Q�C��
G = 70;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�@���ߍ��݁@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���摜�̏���<E1>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�摜�ǂݍ���
ImgNAME = 'BOAT.bmp';
% ImgNAME = 'BARBARA.bmp';
% ImgNAME = 'girl.bmp';
% ImgNAME = 'Airplane.bmp';
% ImgNAME = 'Cameraman.bmp';
% ImgNAME = 'LAX.bmp';
% ImgNAME = 'BRIDGE.bmp';
% ImgNAME = 'Earth.bmp';
% ImgNAME = 'Mandrill.bmp';
% ImgNAME = 'Milkdrop.bmp';
% ImgNAME = 'Pepper.bmp';
% ImgNAME = 'Parrots.bmp';
I = im2double(imread(ImgNAME));

%�\��
figure,imshow(I),title('�z�X�g�摜')

%�����������킹�邽�߂�255�{
I = I * 255;

%�C���[�W���ɂ��� 8 �s 8 ��̃u���b�N�� 2 ���� DCT ���v�Z
T = dctmtx(8);
dct = @(block_struct) T * block_struct.data * T';
dct = blockproc(I,[8 8],dct);

%8x8�̃T�u�u���b�N(1024��)C�ɕ���
C =  mat2cell(dct,[8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8],[8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8]);

%�T�u�u���b�NC�����X�^�X�L�����̏��ŕ��ёւ�
CR = C.';
CR = CR(:)';

%�T�u�u���b�N�����̎擾
Num_block = numel(CR);

%�T�u�u���b�NCR���ƂɃW�O�U�O�X�L����
for i= 1:Num_block
    CRZ{i} = zigzag(CR{i});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�V�[�P���V�̒��o�ƘA��<E2>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�e�T�u�u���b�N��ss�񂩂�se��(����Leg)�𒊏o
se = ss + Leg -1;
for i= 1:Num_block
    seq{i} = CRZ{i}(:,ss:se);
end

%���o�����f�[�^��A��
y = cell2mat(seq);

%���ɖ߂�
I = I / 255;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ߍ��ݏ��̏���<E3>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�t�@�C�����番�U�l���擾
switch Leg
    case 2
        dataHEX_emb = cell2mat((importdata('64.txt')));
    case 4
        dataHEX_emb = cell2mat((importdata('128.txt')));
    case 8
        dataHEX_emb = cell2mat((importdata('256.txt')));
    case 16
        dataHEX_emb = cell2mat((importdata('512.txt')));
    case 32
        dataHEX_emb = cell2mat((importdata('1024.txt')));
    case 50
        dataHEX_emb = cell2mat((importdata('1600.txt')));
    otherwise
        dataHEX_emb = cell2mat((importdata('Emb.txt')));
end

%���ߍ��ރf�[�^���̎擾
Num_data = numel(dataHEX_emb);

%16�i������10�i���֕ϊ�
for i=1:Num_data
    dataDEC(i) = hex2dec(dataHEX_emb(i));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ߍ��ݏ����@<E4> M�L�n��̖��ߍ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M�L�n��(m��)����{����
m = 4;
Md = [mls(m,1),-1];
MdLeg = 2^m;  %M�L�n��̒���

Num_SecInf = 1;        %���������̔ԍ��Bb��(2^m)+1�ɂȂ�����X�V
cnt_emb = 1;        %�n������̃J�E���g(b=2^m�ň���B��������玟�̃V�t�g�����n��̖��ߍ���)

%M�L�n����V�t�g
MdSft = circshift(Md,dataDEC(Num_SecInf));

%���ߍ���
for j=1:Num_data * MdLeg
    if(cnt_emb == MdLeg + 1)
        cnt_emb = 1;
        Num_SecInf = Num_SecInf + 1;
        MdSft = circshift(Md,dataDEC(Num_SecInf));
    end
    Y(j) = y(j) + G * MdSft(cnt_emb); 
    cnt_emb = cnt_emb + 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ߍ��ݏ����A<E5> ����������摜�̐���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%<E2>�Œ��o�����V�[�P���V�ւ̖��ߍ���
j = 1;
iCRZ = CRZ;
for i= 1:Num_block
        iCRZ{i}(:,ss:se) = Y(:,j:j+Leg-1);
        j = j + Leg;
end

%�T�u�u���b�N���Ƃɋt�W�O�U�O�X�L����
for i= 1:Num_block
    iCR{i} = izigzag(iCRZ{i},8,8);
end

%32�~32�ɍĔz��
iC = reshape(iCR,32,32).';

%256�~256�ɘA��
iCl = cell2mat(iC);

%�e�u���b�N�� 2 �����t DCT ���g�p���ăC���[�W���č\��
invdct = @(block_struct) T' * block_struct.data * T;
Iwm = blockproc(iCl,[8 8],invdct);
Iwm = Iwm / 255;

%��������
imwrite(Iwm,'Watermarked_Image_v1.bmp','bmp');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�@���o�@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���o����<D1> �O����<E1><E2>�̎��s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�摜�ǂݍ���
Iwm = im2double(imread('Watermarked_Image_v1.bmp'));

%�����������킹�邽�߂�255�{
Iwm = Iwm * 255;

%�C���[�W���ɂ��� 8 �s 8 ��̃u���b�N�� 2 ���� DCT ���v�Z
dct = @(block_struct) T * block_struct.data * T';
dct = blockproc(Iwm,[8 8],dct);

%8x8�̃T�u�u���b�N(1024��)C�ɕ���
Cx =  mat2cell(dct,[8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8],[8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8]);

%�T�u�u���b�NC�����X�^�X�L�����̏��ŕ��ёւ�
CRx = Cx.';
CRx = CRx(:)';

%�T�u�u���b�NCR���ƂɃW�O�U�O�X�L����
for i= 1:Num_block
    CRZx{i} = zigzag(CRx{i});
end

%�e�T�u�u���b�N��ss�񂩂�se��𒊏o
for i= 1:Num_block
    seq{i} = CRZx{i}(:,ss:se);
end

%���o�����f�[�^��A��
yx = cell2mat(seq);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���o����<D2> ���������̒��o%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�閧���͑��֊֐��̃s�[�N�ʒu�Ƃ��Ă���B
%���ߍ��݂Ɏg����M�L�n��ƒ��o�����f�[�^�̑��ւ����B
j = 1;
for i=1:Num_data
    [acor,lag] = xcorr(yx(:,j:j+MdLeg-1),Md);
    corr{i} = acor;
    j = j + MdLeg;
end

%���ւ��ő�ƂȂ��Ă���ʒu�̒��o
for i=1:Num_data
    [Maximum,Index_M] = max(corr{i});
    if(Index_M >= MdLeg)
        Index_M = Index_M - MdLeg;
    end
    Index(i) = Index_M;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���o����<D3> �č\��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%10�i���̈ʒu�f�[�^��16�i���֕ϊ���string�^�œ����������č\��
for i=1:Num_data
    dataHEX_ext(i) = lower(string(dec2hex(Index(i))));
end

%���ɖ߂�
Iwm = Iwm / 255;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  �]��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSNR
PSNR = psnr(Iwm,I);
fprintf('\n PSNR : %0.2f[dB]\n', PSNR);

%SSIM
SSIM = ssim(Iwm,I);
fprintf('\n SSIM : %0.4f\n', SSIM);

%BER
[number,ratio] = biterr(Index,dataDEC);
fprintf('\n BER  : %0.4f[��]\n', ratio*100);

%����������摜�̕\��
figure,imshow(Iwm),title('����������摜')

%�f�[�^���e�L�X�g�t�@�C���֏o��
txt = fopen('Ext.txt','w');
fprintf(txt,'Image : %s',ImgNAME);
fprintf(txt,'\n');
fprintf(txt,'Embedding start number : %d\n',ss);
fprintf(txt,'Gain : %0.3f\n',G);
fprintf(txt,'Used M�Lsequence : ');
fprintf(txt,'%d',Md);
fprintf(txt,'\n\n');
fprintf(txt,'Embedding_data   : ');
fprintf(txt,'%s',dataHEX_emb);
fprintf(txt,'\n');
fprintf(txt,'Extraction_data  : ');
fprintf(txt,'%s',dataHEX_ext);
fprintf(txt,'\n\n');
fprintf(txt,'PSNR : %0.2f[dB]\n',PSNR');
fprintf(txt,'SSIM : %0.2f[dB]\n',SSIM');
fprintf(txt,'BER  : %d[bits] %0.4f[��]\n',number,ratio*100);
fclose(txt);
fprintf('\n�f�[�^��Ext.txt�ɏo�͂��܂����B\n');


