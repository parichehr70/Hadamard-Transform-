clear all; close all; clc;
I = im2double(imread('cameraman.tif'));
N = 256;  % Length of Walsh (Hadamard) functions
hadamardMatrix = hadamard(N);
HadIdx = 0:N-1;                          % Hadamard index
M = log2(N)+1;                           % Number of bits to represent the index
binHadIdx = fliplr(dec2bin(HadIdx,M))-'0'; % Bit reversing of the binary index
binSeqIdx = zeros(N,M-1);                  % Pre-allocate memory
for k = M:-1:2
    % Binary sequency index
    binSeqIdx(:,k) = xor(binHadIdx(:,k),binHadIdx(:,k-1));
end
SeqIdx = binSeqIdx*pow2((M-1:-1:0)');    % Binary to integer sequency index
W = hadamardMatrix(SeqIdx+1,:); % 1-based indexing
J = W*I*W;
n = 512;
[e,f] = meshgrid(-1:1/((n-1)/2):1,-1:1/((n-1)/2):1);
LOW1 = 1/(2*pi*10^2)*exp(-(e.^2+f.^2)/2*10^2);       % sigma=10
figure,mesh(LOW1);
LOW1 = LOW1(257:512,257:512);
figure,mesh(LOW1);
LOW1 = (LOW1-min(LOW1(:)))/(max(LOW1(:)-min(LOW1(:))));
LOW2 = 1/(2*pi*5^2)*exp(-(e.^2+f.^2)/2*5^2);       % sigma=5
LOW2 = LOW2(257:512,257:512);
LOW2 = (LOW2-min(LOW2(:)))/(max(LOW2(:)-min(LOW2(:))));
LOW3 = 1/(2*pi*0.1^2)*exp(-(e.^2+f.^2)/2*0.1^2);      % sigma=0.1
LOW3 = LOW3(257:512,257:512);
LOW3 = (LOW3-min(LOW3(:)))/(max(LOW3(:)-min(LOW3(:))));
HI1 = 1 - LOW1;HI1 = (HI1-min(HI1(:)))/(max(HI1(:)-min(HI1(:))));
HI2 = 1 - LOW2;HI2 = (HI2-min(HI2(:)))/(max(HI2(:)-min(HI2(:))));
HI3 = 1 - LOW3;HI3 = (HI3-min(HI3(:)))/(max(HI3(:)-min(HI3(:))));
BS1 = LOW1 + HI2 ; BS1 = (BS1-min(BS1(:)))/(max(BS1(:)-min(BS1(:))));
BS2 = LOW2 + HI3 ; BS2 = (BS2-min(BS2(:)))/(max(BS2(:)-min(BS2(:))));
BP1 = 1 - BS1 ; BP1 = (BP1-min(BP1(:)))/(max(BP1(:)-min(BP1(:))));
BP2 = 1 - BS2 ; BP2 = (BP2-min(BP2(:)))/(max(BP2(:)-min(BP2(:))));
YL1 = LOW1.*J ; figure,imshow(W*YL1*W/256^2);title('LowPass 1');
YL2 = LOW2.*J ; figure,imshow(W*YL2*W/256^2);title('LowPass 2');
YH1 = HI1.*J ; figure,imshow(W*YH1*W/256^2);title('HighPass 1');
YH2 = HI2.*J ; figure,imshow(W*YH2*W/256^2);title('HighPass 2');
YBP1 = BP1.*J ; figure,imshow(W*YBP1*W/256^2);title('BandPass 1');
YBP2 = BP2.*J ; figure,imshow(W*YBP2*W/256^2);title('BandPass 2');
YBS1 = BS1.*J ; figure,imshow(W*YBS1*W/256^2);title('BandStop 1');
YBS2 = BS2.*J ; figure,imshow(W*YBS2*W/256^2);title('BandStop 2');
