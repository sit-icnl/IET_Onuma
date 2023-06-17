
%��ҁF����C�m�@'20�� af16023
%               '22�� ma20017

% ver2 RS������p���Č��������������@



clc;
clear all;
close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�@��v�p�����[�^�ݒ�@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ߍ��݊J�n�ʒu
ss = 32;  %MAX64-Leg+1

%RS����RS(R,S) ������R���w�肷��΂��ׂĎ����I�Ƀp���[���^�����߂�悤�ɐݒ�ς�
%R = 15;
%R = 31;
%R = 63;
%R = 127;
R = 255;
%R = 511;

%1�u���b�N�������DCT�W���̒��o��(���ɉ����Ė��ߍ��ޏ��ʂ��ő���܂œ���Ȃ��ƃG���[�D���������ς�)
%Leg = 2; %64
%Leg = 4; %128
%Leg = 8; %256
Leg = 16; %512
%Leg = 32; %1024
%Leg = 64; %2048

%M'�n��̃Q�C��
G = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�@���ߍ��݁@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���摜�̏���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�V�[�P���V�̒��o�ƘA��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�e�T�u�u���b�N��ss�񂩂�se��(����Leg)�𒊏o
se = ss + Leg -1;
for i= 1:Num_block
    seq{i} = CRZ{i}(:,ss:se);
end

%���o�����f�[�^��A��
y = cell2mat(seq);

%���ɖ߂�
I = I / 255;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ߍ��ݏ��̏���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    case 64
        dataHEX_emb = cell2mat((importdata('2048.txt')));
    otherwise
        dataHEX_emb = cell2mat((importdata('Emb.txt')));
end

%���ߍ��ރf�[�^���̎擾
Num_data = numel(dataHEX_emb);

%16�i������10�i���֕ϊ�
for i=1:Num_data
    dataDEC(i) = hex2dec(dataHEX_emb(i));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ߍ��ݏ����@ M�L�n��̖��ߍ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M�L�n��(m��)����{����
m_s = 4;
Md = [mls(m_s,1),-1];
MdLeg = 2^m_s;  %M�L�n��̒���

Num_SecInfo = 1;        %���������̔ԍ��Bb��(2^m_s)+1�ɂȂ�����X�V
cnt_emb = 1;        %�n������̃J�E���g(b=2^m_s�ň���B��������玟�̃V�t�g�����n��̖��ߍ���)

%M�L�n����V�t�g
MdSft = circshift(Md,dataDEC(Num_SecInfo));

%���ߍ���
for j=1:Num_data * MdLeg
    if(cnt_emb == MdLeg + 1)
        cnt_emb = 1;
        Num_SecInfo = Num_SecInfo + 1;
        MdSft = circshift(Md,dataDEC(Num_SecInfo));
    end
    Y(j) = y(j) + G * MdSft(cnt_emb); 
    cnt_emb = cnt_emb + 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ߍ��ݏ����A ����������摜�̐���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
imwrite(Iwm,'Watermarked_Image_0.bmp','bmp');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RS�����̖��ߍ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ix = imread('Watermarked_Image_0.bmp');

%RS�����Ɋւ���p�����[�^
m = log2(R+1);            % 1�V���{��������̃r�b�g��
n = 2^m-1;                % ������
k = 2^(m-1)-1;            % ���V���{����
q = fix(Num_data/k);      %RS�����̖{�����Z�o(q+1�ŏI��)
i=0; Num_SecInfo=1; cnt_emb=k+1; count=0;

while(1)
    %k�����������Ă����C�]�肪�ł�ꏊ(q+1���)��0�Ńp�f�B���O�D
    count = count + 1;
    if(count <= q)
        for p=1:k
            d(p) = dataDEC(1,i+p);
        end
    end
    if(count == q +1)
        for p=1:Num_data - q*k
            d(p) = dataDEC(1,i+p);
        end
        for p=Num_data - q*k+1:k
            d(p) = 0;
        end
    end
    
    %10�i���̔閧��񂩂�K���A�̔z��̐���(������̐���)
    switch n
    case 15
        msg = gf([d(1) d(2) d(3) d(4) d(5) d(6) d(7)],m);
    case 31
        msg = gf([d(1) d(2) d(3) d(4) d(5) d(6) d(7) d(8) d(9) d(10) d(11) d(12) d(13) d(14) d(15)],m);
    case 63
        msg = gf([d(1) d(2) d(3) d(4) d(5) d(6) d(7) d(8) d(9) d(10) d(11) d(12) d(13) d(14) d(15) d(16) d(17) d(18) d(19) d(20) d(21) d(22) d(23) d(24) d(25) d(26) d(27) d(28) d(29) d(30) d(31)],m);
    case 127
        msg = gf([d(1) d(2) d(3) d(4) d(5) d(6) d(7) d(8) d(9) d(10) d(11) d(12) d(13) d(14) d(15) d(16) d(17) d(18) d(19) d(20) d(21) d(22) d(23) d(24) d(25) d(26) d(27) d(28) d(29) d(30) d(31) d(32) d(33) d(34) d(35) d(36) d(37) d(38) d(39) d(40) d(41) d(42) d(43) d(44) d(45) d(46) d(47) d(48) d(49) d(50) d(51) d(52) d(53) d(54) d(55) d(56) d(57) d(58) d(59) d(60) d(61) d(62) d(63)],m);
    case 255
        msg = gf([d(1) d(2) d(3) d(4) d(5) d(6) d(7) d(8) d(9) d(10) d(11) d(12) d(13) d(14) d(15) d(16) d(17) d(18) d(19) d(20) d(21) d(22) d(23) d(24) d(25) d(26) d(27) d(28) d(29) d(30) d(31) d(32) d(33) d(34) d(35) d(36) d(37) d(38) d(39) d(40) d(41) d(42) d(43) d(44) d(45) d(46) d(47) d(48) d(49) d(50) d(51) d(52) d(53) d(54) d(55) d(56) d(57) d(58) d(59) d(60) d(61) d(62) d(63) d(64) d(65) d(66) d(67) d(68) d(69) d(70) d(71) d(72) d(73) d(74) d(75) d(76) d(77) d(78) d(79) d(80) d(81) d(82) d(83) d(84) d(85) d(86) d(87) d(88) d(89) d(90) d(91) d(92) d(93) d(94) d(95) d(96) d(97) d(98) d(99) d(100) d(101) d(102) d(103) d(104) d(105) d(106) d(107) d(108) d(109) d(110) d(111) d(112) d(113) d(114) d(115) d(116) d(117) d(118) d(119) d(120) d(121) d(122) d(123) d(124) d(125) d(126) d(127)],m);
    case 511
        msg = gf([d(1) d(2) d(3) d(4) d(5) d(6) d(7) d(8) d(9) d(10) d(11) d(12) d(13) d(14) d(15) d(16) d(17) d(18) d(19) d(20) d(21) d(22) d(23) d(24) d(25) d(26) d(27) d(28) d(29) d(30) d(31) d(32) d(33) d(34) d(35) d(36) d(37) d(38) d(39) d(40) d(41) d(42) d(43) d(44) d(45) d(46) d(47) d(48) d(49) d(50) d(51) d(52) d(53) d(54) d(55) d(56) d(57) d(58) d(59) d(60) d(61) d(62) d(63) d(64) d(65) d(66) d(67) d(68) d(69) d(70) d(71) d(72) d(73) d(74) d(75) d(76) d(77) d(78) d(79) d(80) d(81) d(82) d(83) d(84) d(85) d(86) d(87) d(88) d(89) d(90) d(91) d(92) d(93) d(94) d(95) d(96) d(97) d(98) d(99) d(100) d(101) d(102) d(103) d(104) d(105) d(106) d(107) d(108) d(109) d(110) d(111) d(112) d(113) d(114) d(115) d(116) d(117) d(118) d(119) d(120) d(121) d(122) d(123) d(124) d(125) d(126) d(127) d(128) d(129) d(130) d(131) d(132) d(133) d(134) d(135) d(136) d(137) d(138) d(139) d(140) d(141) d(142) d(143) d(144) d(145) d(146) d(147) d(148) d(149) d(150) d(151) d(152) d(153) d(154) d(155) d(156) d(157) d(158) d(159) d(160) d(161) d(162) d(163) d(164) d(165) d(166) d(167) d(168) d(169) d(170) d(171) d(172) d(173) d(174) d(175) d(176) d(177) d(178) d(179) d(180) d(181) d(182) d(183) d(184) d(185) d(186) d(187) d(188) d(189) d(190) d(191) d(192) d(193) d(194) d(195) d(196) d(197) d(198) d(199) d(200) d(201) d(202) d(203) d(204) d(205) d(206) d(207) d(208) d(209) d(210) d(211) d(212) d(213) d(214) d(215) d(216) d(217) d(218) d(219) d(220) d(221) d(222) d(223) d(224) d(225) d(226) d(227) d(228) d(229) d(230) d(231) d(232) d(233) d(234) d(235) d(236) d(237) d(238) d(239) d(240) d(241) d(242) d(243) d(244) d(245) d(246) d(247) d(248) d(249) d(250) d(251) d(252) d(253) d(254) d(255)],m);
    otherwise
        beep
        fprintf('�G���[�FR���s���Ȓl�D\n');
        return
    end
   
    code = rsenc(msg,n,k);
    
    %�����ꂩ��璷�V���{���̒��o
    Ary(:,Num_SecInfo:cnt_emb) = code.x(:,k+1:n);
    
    i = i + k;
    Num_SecInfo = Num_SecInfo + k+1;
    cnt_emb = cnt_emb + k+1;
    
    if(count == q + 1)
        break;
    end
end
%����������2�i����
binAry = de2bi(Ary);

%����������(���ߍ��ސ�)�̎擾
numAry = numel(binAry);

%�z��̃T�C�Y���擾
sz = size(binAry);

%���ɂ܂Ƃ�
binAry = binAry.';
binAry = binAry(:);

%��f�l��2�i���֕ϊ�
binStr = de2bi(Ix);
binStr2 = de2bi(Ix);

%�摜�̍ŉ��ʃr�b�g�֒��������𖄂ߍ���
for i=1:numAry
    binStr2(i,1) = binAry(i,1);
end

%2�i������10�i���֖߂�
decStr = bi2de(binStr2);

%�摜�̍č\�z
z=1; u=1; h=256;
while(1)
    Ix(1:256,z) = decStr(u:h);
    u=u+256; h=h+256;
if(z == 256)
    break;
end
z = z+1;
end
imwrite(Ix,'Watermarked_Image_v2.bmp','bmp');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�@���o�@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���o���� �O�����̎��s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�摜�ǂݍ���
Ix1 = im2double(imread('Watermarked_Image_v2.bmp'));
Ix2 = imread('Watermarked_Image_v2.bmp');

%�����������킹�邽�߂�255�{
Ix1 = Ix1 * 255;

%�C���[�W���ɂ��� 8 �s 8 ��̃u���b�N�� 2 ���� DCT ���v�Z
dct = @(block_struct) T * block_struct.data * T';
dct = blockproc(Ix1,[8 8],dct);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���o���� ���������̒��o%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%�@������̒��o %%%%%%%%%%%%%%%%%%%
%��f�l��2�i���֕ϊ�
binStrX = de2bi(Ix2);
codeX = binStrX(1:numAry,1);
codeX = reshape(codeX,sz(1,2),sz(1,1)).';
codeX = bi2de(codeX).';

count=0; u=1; i=0; j=0; w=1;

%�摜���璊�o�����閧��� In �ƌ��������� dx �̓���
while(1)
    count = count + 1;
    if(count <= q)
        for p=1:k+1
            if(p<=k)
                In(p) = Index(1,j+p);
                dx(p) = codeX(1,i+p);
            end
            if(p>k)
                dx(p) = codeX(1,i+p);
            end
        end
    end
    
    if(count == q + 1)
        for p=1:Num_data - q*k
            In(p) = Index(1,j+p);
            dx(p) = codeX(1,i+p);
        end
        for p=Num_data - q*k+1:k+1
            In(p) = 0;
            dx(p) = codeX(1,i+p);
        end
    end

    %In��dx�𓝍������K���A�̔z��̐���
    switch n
    case 15
        msgX = gf([In(1) In(2) In(3) In(4) In(5) In(6) In(7) dx(1) dx(2) dx(3) dx(4) dx(5) dx(6) dx(7) dx(8)],m);
    case 31
        msgX = gf([In(1) In(2) In(3) In(4) In(5) In(6) In(7) In(8) In(9) In(10) In(11) In(12) In(13) In(14) In(15) dx(1) dx(2) dx(3) dx(4) dx(5) dx(6) dx(7) dx(8) dx(9) dx(10) dx(11) dx(12) dx(13) dx(14) dx(15) dx(16)],m);
    case 63
        msgX = gf([In(1) In(2) In(3) In(4) In(5) In(6) In(7) In(8) In(9) In(10) In(11) In(12) In(13) In(14) In(15) In(16) In(17) In(18) In(19) In(20) In(21) In(22) In(23) In(24) In(25) In(26) In(27) In(28) In(29) In(30) In(31) dx(1) dx(2) dx(3) dx(4) dx(5) dx(6) dx(7) dx(8) dx(9) dx(10) dx(11) dx(12) dx(13) dx(14) dx(15) dx(16) dx(17) dx(18) dx(19) dx(20) dx(21) dx(22) dx(23) dx(24) dx(25) dx(26) dx(27) dx(28) dx(29) dx(30) dx(31) dx(32)],m);
    case 127
        msgX = gf([In(1) In(2) In(3) In(4) In(5) In(6) In(7) In(8) In(9) In(10) In(11) In(12) In(13) In(14) In(15) In(16) In(17) In(18) In(19) In(20) In(21) In(22) In(23) In(24) In(25) In(26) In(27) In(28) In(29) In(30) In(31) In(32) In(33) In(34) In(35) In(36) In(37) In(38) In(39) In(40) In(41) In(42) In(43) In(44) In(45) In(46) In(47) In(48) In(49) In(50) In(51) In(52) In(53) In(54) In(55) In(56) In(57) In(58) In(59) In(60) In(61) In(62) In(63) dx(1) dx(2) dx(3) dx(4) dx(5) dx(6) dx(7) dx(8) dx(9) dx(10) dx(11) dx(12) dx(13) dx(14) dx(15) dx(16) dx(17) dx(18) dx(19) dx(20) dx(21) dx(22) dx(23) dx(24) dx(25) dx(26) dx(27) dx(28) dx(29) dx(30) dx(31) dx(32) dx(33) dx(34) dx(35) dx(36) dx(37) dx(38) dx(39) dx(40) dx(41) dx(42) dx(43) dx(44) dx(45) dx(46) dx(47) dx(48) dx(49) dx(50) dx(51) dx(52) dx(53) dx(54) dx(55) dx(56) dx(57) dx(58) dx(59) dx(60) dx(61) dx(62) dx(63) dx(64)],m);
    case 255
        msgX = gf([In(1) In(2) In(3) In(4) In(5) In(6) In(7) In(8) In(9) In(10) In(11) In(12) In(13) In(14) In(15) In(16) In(17) In(18) In(19) In(20) In(21) In(22) In(23) In(24) In(25) In(26) In(27) In(28) In(29) In(30) In(31) In(32) In(33) In(34) In(35) In(36) In(37) In(38) In(39) In(40) In(41) In(42) In(43) In(44) In(45) In(46) In(47) In(48) In(49) In(50) In(51) In(52) In(53) In(54) In(55) In(56) In(57) In(58) In(59) In(60) In(61) In(62) In(63) In(64) In(65) In(66) In(67) In(68) In(69) In(70) In(71) In(72) In(73) In(74) In(75) In(76) In(77) In(78) In(79) In(80) In(81) In(82) In(83) In(84) In(85) In(86) In(87) In(88) In(89) In(90) In(91) In(92) In(93) In(94) In(95) In(96) In(97) In(98) In(99) In(100) In(101) In(102) In(103) In(104) In(105) In(106) In(107) In(108) In(109) In(110) In(111) In(112) In(113) In(114) In(115) In(116) In(117) In(118) In(119) In(120) In(121) In(122) In(123) In(124) In(125) In(126) In(127) dx(1) dx(2) dx(3) dx(4) dx(5) dx(6) dx(7) dx(8) dx(9) dx(10) dx(11) dx(12) dx(13) dx(14) dx(15) dx(16) dx(17) dx(18) dx(19) dx(20) dx(21) dx(22) dx(23) dx(24) dx(25) dx(26) dx(27) dx(28) dx(29) dx(30) dx(31) dx(32) dx(33) dx(34) dx(35) dx(36) dx(37) dx(38) dx(39) dx(40) dx(41) dx(42) dx(43) dx(44) dx(45) dx(46) dx(47) dx(48) dx(49) dx(50) dx(51) dx(52) dx(53) dx(54) dx(55) dx(56) dx(57) dx(58) dx(59) dx(60) dx(61) dx(62) dx(63) dx(64) dx(65) dx(66) dx(67) dx(68) dx(69) dx(70) dx(71) dx(72) dx(73) dx(74) dx(75) dx(76) dx(77) dx(78) dx(79) dx(80) dx(81) dx(82) dx(83) dx(84) dx(85) dx(86) dx(87) dx(88) dx(89) dx(90) dx(91) dx(92) dx(93) dx(94) dx(95) dx(96) dx(97) dx(98) dx(99) dx(100) dx(101) dx(102) dx(103) dx(104) dx(105) dx(106) dx(107) dx(108) dx(109) dx(110) dx(111) dx(112) dx(113) dx(114) dx(115) dx(116) dx(117) dx(118) dx(119) dx(120) dx(121) dx(122) dx(123) dx(124) dx(125) dx(126) dx(127) dx(128)],m);
    case 511
        msgX = gf([In(1) In(2) In(3) In(4) In(5) In(6) In(7) In(8) In(9) In(10) In(11) In(12) In(13) In(14) In(15) In(16) In(17) In(18) In(19) In(20) In(21) In(22) In(23) In(24) In(25) In(26) In(27) In(28) In(29) In(30) In(31) In(32) In(33) In(34) In(35) In(36) In(37) In(38) In(39) In(40) In(41) In(42) In(43) In(44) In(45) In(46) In(47) In(48) In(49) In(50) In(51) In(52) In(53) In(54) In(55) In(56) In(57) In(58) In(59) In(60) In(61) In(62) In(63) In(64) In(65) In(66) In(67) In(68) In(69) In(70) In(71) In(72) In(73) In(74) In(75) In(76) In(77) In(78) In(79) In(80) In(81) In(82) In(83) In(84) In(85) In(86) In(87) In(88) In(89) In(90) In(91) In(92) In(93) In(94) In(95) In(96) In(97) In(98) In(99) In(100) In(101) In(102) In(103) In(104) In(105) In(106) In(107) In(108) In(109) In(110) In(111) In(112) In(113) In(114) In(115) In(116) In(117) In(118) In(119) In(120) In(121) In(122) In(123) In(124) In(125) In(126) In(127) In(128) In(129) In(130) In(131) In(132) In(133) In(134) In(135) In(136) In(137) In(138) In(139) In(140) In(141) In(142) In(143) In(144) In(145) In(146) In(147) In(148) In(149) In(150) In(151) In(152) In(153) In(154) In(155) In(156) In(157) In(158) In(159) In(160) In(161) In(162) In(163) In(164) In(165) In(166) In(167) In(168) In(169) In(170) In(171) In(172) In(173) In(174) In(175) In(176) In(177) In(178) In(179) In(180) In(181) In(182) In(183) In(184) In(185) In(186) In(187) In(188) In(189) In(190) In(191) In(192) In(193) In(194) In(195) In(196) In(197) In(198) In(199) In(200) In(201) In(202) In(203) In(204) In(205) In(206) In(207) In(208) In(209) In(210) In(211) In(212) In(213) In(214) In(215) In(216) In(217) In(218) In(219) In(220) In(221) In(222) In(223) In(224) In(225) In(226) In(227) In(228) In(229) In(230) In(231) In(232) In(233) In(234) In(235) In(236) In(237) In(238) In(239) In(240) In(241) In(242) In(243) In(244) In(245) In(246) In(247) In(248) In(249) In(250) In(251) In(252) In(253) In(254) In(255) dx(1) dx(2) dx(3) dx(4) dx(5) dx(6) dx(7) dx(8) dx(9) dx(10) dx(11) dx(12) dx(13) dx(14) dx(15) dx(16) dx(17) dx(18) dx(19) dx(20) dx(21) dx(22) dx(23) dx(24) dx(25) dx(26) dx(27) dx(28) dx(29) dx(30) dx(31) dx(32) dx(33) dx(34) dx(35) dx(36) dx(37) dx(38) dx(39) dx(40) dx(41) dx(42) dx(43) dx(44) dx(45) dx(46) dx(47) dx(48) dx(49) dx(50) dx(51) dx(52) dx(53) dx(54) dx(55) dx(56) dx(57) dx(58) dx(59) dx(60) dx(61) dx(62) dx(63) dx(64) dx(65) dx(66) dx(67) dx(68) dx(69) dx(70) dx(71) dx(72) dx(73) dx(74) dx(75) dx(76) dx(77) dx(78) dx(79) dx(80) dx(81) dx(82) dx(83) dx(84) dx(85) dx(86) dx(87) dx(88) dx(89) dx(90) dx(91) dx(92) dx(93) dx(94) dx(95) dx(96) dx(97) dx(98) dx(99) dx(100) dx(101) dx(102) dx(103) dx(104) dx(105) dx(106) dx(107) dx(108) dx(109) dx(110) dx(111) dx(112) dx(113) dx(114) dx(115) dx(116) dx(117) dx(118) dx(119) dx(120) dx(121) dx(122) dx(123) dx(124) dx(125) dx(126) dx(127) dx(128) dx(129) dx(130) dx(131) dx(132) dx(133) dx(134) dx(135) dx(136) dx(137) dx(138) dx(139) dx(140) dx(141) dx(142) dx(143) dx(144) dx(145) dx(146) dx(147) dx(148) dx(149) dx(150) dx(151) dx(152) dx(153) dx(154) dx(155) dx(156) dx(157) dx(158) dx(159) dx(160) dx(161) dx(162) dx(163) dx(164) dx(165) dx(166) dx(167) dx(168) dx(169) dx(170) dx(171) dx(172) dx(173) dx(174) dx(175) dx(176) dx(177) dx(178) dx(179) dx(180) dx(181) dx(182) dx(183) dx(184) dx(185) dx(186) dx(187) dx(188) dx(189) dx(190) dx(191) dx(192) dx(193) dx(194) dx(195) dx(196) dx(197) dx(198) dx(199) dx(200) dx(201) dx(202) dx(203) dx(204) dx(205) dx(206) dx(207) dx(208) dx(209) dx(210) dx(211) dx(212) dx(213) dx(214) dx(215) dx(216) dx(217) dx(218) dx(219) dx(220) dx(221) dx(222) dx(223) dx(224) dx(225) dx(226) dx(227) dx(228) dx(229) dx(230) dx(231) dx(232) dx(233) dx(234) dx(235) dx(236) dx(237) dx(238) dx(239) dx(240) dx(241) dx(242) dx(243) dx(244) dx(245) dx(246) dx(247) dx(248) dx(249) dx(250) dx(251) dx(252) dx(253) dx(254) dx(255) dx(256)],m);
    otherwise
        return
    end
        
    %����
    [rxcode,cnumerr,code] = rsdec(msgX,n,k);
    fprintf('\n Number of corrected information  [%d] : %d\n',w,cnumerr);
    
    %�����������U�l���i�[
    fixed(1,u:u+k-1) = rxcode.x(1,1:k);
    
    i=i+k+1; u=u+k; j=j+k; w=w+1;
    if(count == q +1)
        break;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���o���� �č\��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%10�i���̈ʒu�f�[�^��16�i���֕ϊ���string�^�œ����������č\��
for i=1:Num_data
    dataHEX_ext(i) = lower(string(dec2hex(fixed(i))));
end

%���ɖ߂�
Ix1 = Ix1 / 255;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  �]��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSNR
PSNR = psnr(Ix1,I);
fprintf('\n PSNR : %0.2f[dB]\n', PSNR);

%SSIM
SSIM = ssim(Ix1,I);
fprintf('\n SSIM : %0.4f\n', SSIM);

%BER�@�������O
[number1,ratio1] = biterr(Index,dataDEC);      
fprintf('\n BER (Before correction) : %0.4f[��]\n', ratio1*100);

%BER�@������
fixed(:,Num_data+1:(q+1)*k) = [];  %0�Ńp�f�B���O���������폜
fixed = double(fixed);
[number2,ratio2] = biterr(fixed,dataDEC);
fprintf('\n BER (After correction)  : %0.4f[��]\n', ratio2*100);

%����������摜�̕\��
figure,imshow(Ix1),title('����������摜')

%�f�[�^���e�L�X�g�t�@�C���֏o��
% txt = fopen('Ext.txt','w');
% fprintf(txt,'Image : %s',ImgNAME);
% fprintf(txt,'\n');
% fprintf(txt,'Embedding start number : %d\n',ss);
% fprintf(txt,'Gain : %0.3f\n',G);
% fprintf(txt,'Used M�Lsequence : ');
% fprintf(txt,'%d',Md);
% fprintf(txt,'\n\n');
% fprintf(txt,'Embedding_data                      : ');
% fprintf(txt,'%s',dataHEX_emb);
% fprintf(txt,'\n');
% fprintf(txt,'Extraction_data (Before correction) : ');
% fprintf(txt,'%s',lower(string(dec2hex(Index))));
% fprintf(txt,'\n');
% fprintf(txt,'Extraction_data (After correction)  : ');
% fprintf(txt,'%s',dataHEX_ext);
% fprintf(txt,'\n\n');
% fprintf(txt,'PSNR : %0.2f[dB]\n',PSNR');
% fprintf(txt,'SSIM : %0.2f[dB]\n',SSIM');
% fprintf(txt,'BER (Before correction) : %d[bits] %0.4f[��]\n',number1,ratio1*100);
% fprintf(txt,'BER (After correction)  : %d[bits] %0.4f[��]\n',number2,ratio2*100);
% fclose(txt);
% fprintf('\n�f�[�^��Ext.txt�ɏo�͂��܂����B\n');


