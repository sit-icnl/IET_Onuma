
%作者：大沼海仁　'20卒 af16023
%               '22卒 ma20017

%スペクトル拡散法と(k,n)閾値秘密分散法を用いたステガノグラフィ
%暗号化サイト：http://point-at-infinity.org/ssss/demo.html




clc;
clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　主要パラメータ設定　%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%埋め込み開始位置
ss = 5;  %MAX64-Leg+1

%DCT係数の抽出数(数に応じて埋め込む情報量を最大限まで入れないとエラー)
Leg = 2; %64
%Leg = 4; %128
%Leg = 8; %256
%Leg = 16; %512
%Leg = 32; %1024
%Leg = 50; %1600

%M'系列のゲイン
%G = 4.159;   %ゲイン
%G = 7.699;   %ゲイン
%G = 11.036;  %ゲイン
%G = 19.127;  %ゲイン
G = 70;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　埋め込み　%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%原画像の処理<E1>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画像読み込み
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

%表示
figure,imshow(I),title('ホスト画像')

%整合性を合わせるために255倍
I = I * 255;

%イメージ内にある 8 行 8 列のブロックの 2 次元 DCT を計算
T = dctmtx(8);
dct = @(block_struct) T * block_struct.data * T';
dct = blockproc(I,[8 8],dct);

%8x8のサブブロック(1024個)Cに分割
C =  mat2cell(dct,[8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8],[8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8]);

%サブブロックCをラスタスキャンの順で並び替え
CR = C.';
CR = CR(:)';

%サブブロック総数の取得
Num_block = numel(CR);

%サブブロックCRごとにジグザグスキャン
for i= 1:Num_block
    CRZ{i} = zigzag(CR{i});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%シーケンシの抽出と連結<E2>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%各サブブロックのss列からse列(長さLeg)を抽出
se = ss + Leg -1;
for i= 1:Num_block
    seq{i} = CRZ{i}(:,ss:se);
end

%抽出したデータを連結
y = cell2mat(seq);

%元に戻す
I = I / 255;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%埋め込み情報の処理<E3>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ファイルから分散値を取得
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

%埋め込むデータ数の取得
Num_data = numel(dataHEX_emb);

%16進数から10進数へ変換
for i=1:Num_data
    dataDEC(i) = hex2dec(dataHEX_emb(i));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%埋め込み処理�@<E4> M´系列の埋め込み%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M´系列(m次)を一本生成
m = 4;
Md = [mls(m,1),-1];
MdLeg = 2^m;  %M´系列の長さ

Num_SecInf = 1;        %透かし情報の番号。bが(2^m)+1になったら更新
cnt_emb = 1;        %系列周期のカウント(b=2^mで一周。一周したら次のシフトした系列の埋め込み)

%M´系列をシフト
MdSft = circshift(Md,dataDEC(Num_SecInf));

%埋め込み
for j=1:Num_data * MdLeg
    if(cnt_emb == MdLeg + 1)
        cnt_emb = 1;
        Num_SecInf = Num_SecInf + 1;
        MdSft = circshift(Md,dataDEC(Num_SecInf));
    end
    Y(j) = y(j) + G * MdSft(cnt_emb); 
    cnt_emb = cnt_emb + 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%埋め込み処理�A<E5> 透かし入り画像の生成%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%<E2>で抽出したシーケンシへの埋め込み
j = 1;
iCRZ = CRZ;
for i= 1:Num_block
        iCRZ{i}(:,ss:se) = Y(:,j:j+Leg-1);
        j = j + Leg;
end

%サブブロックごとに逆ジグザグスキャン
for i= 1:Num_block
    iCR{i} = izigzag(iCRZ{i},8,8);
end

%32×32に再配列
iC = reshape(iCR,32,32).';

%256×256に連結
iCl = cell2mat(iC);

%各ブロックの 2 次元逆 DCT を使用してイメージを再構成
invdct = @(block_struct) T' * block_struct.data * T;
Iwm = blockproc(iCl,[8 8],invdct);
Iwm = Iwm / 255;

%書き込み
imwrite(Iwm,'Watermarked_Image_v1.bmp','bmp');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%　抽出　%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%抽出処理<D1> 前処理<E1><E2>の実行%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画像読み込み
Iwm = im2double(imread('Watermarked_Image_v1.bmp'));

%整合性を合わせるために255倍
Iwm = Iwm * 255;

%イメージ内にある 8 行 8 列のブロックの 2 次元 DCT を計算
dct = @(block_struct) T * block_struct.data * T';
dct = blockproc(Iwm,[8 8],dct);

%8x8のサブブロック(1024個)Cに分割
Cx =  mat2cell(dct,[8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8],[8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8]);

%サブブロックCをラスタスキャンの順で並び替え
CRx = Cx.';
CRx = CRx(:)';

%サブブロックCRごとにジグザグスキャン
for i= 1:Num_block
    CRZx{i} = zigzag(CRx{i});
end

%各サブブロックのss列からse列を抽出
for i= 1:Num_block
    seq{i} = CRZx{i}(:,ss:se);
end

%抽出したデータを連結
yx = cell2mat(seq);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%抽出処理<D2> 透かし情報の抽出%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%秘密情報は相関関数のピーク位置としている。
%埋め込みに使ったM´系列と抽出したデータの相関を取る。
j = 1;
for i=1:Num_data
    [acor,lag] = xcorr(yx(:,j:j+MdLeg-1),Md);
    corr{i} = acor;
    j = j + MdLeg;
end

%相関が最大となっている位置の抽出
for i=1:Num_data
    [Maximum,Index_M] = max(corr{i});
    if(Index_M >= MdLeg)
        Index_M = Index_M - MdLeg;
    end
    Index(i) = Index_M;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%抽出処理<D3> 再構成%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%10進数の位置データを16進数へ変換しstring型で透かし情報を再構成
for i=1:Num_data
    dataHEX_ext(i) = lower(string(dec2hex(Index(i))));
end

%元に戻す
Iwm = Iwm / 255;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  評価  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSNR
PSNR = psnr(Iwm,I);
fprintf('\n PSNR : %0.2f[dB]\n', PSNR);

%SSIM
SSIM = ssim(Iwm,I);
fprintf('\n SSIM : %0.4f\n', SSIM);

%BER
[number,ratio] = biterr(Index,dataDEC);
fprintf('\n BER  : %0.4f[％]\n', ratio*100);

%透かし入り画像の表示
figure,imshow(Iwm),title('透かし入り画像')

%データをテキストファイルへ出力
txt = fopen('Ext.txt','w');
fprintf(txt,'Image : %s',ImgNAME);
fprintf(txt,'\n');
fprintf(txt,'Embedding start number : %d\n',ss);
fprintf(txt,'Gain : %0.3f\n',G);
fprintf(txt,'Used M´sequence : ');
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
fprintf(txt,'BER  : %d[bits] %0.4f[％]\n',number,ratio*100);
fclose(txt);
fprintf('\nデータをExt.txtに出力しました。\n');


