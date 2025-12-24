% main.m
clc; clear;

sigma_E = 10^5;      % 加密查找表(ELUT)的标准差
T = 3000;            % 加密查找表(ELUT)的长度
R = 4;               % 每个投影方向对应的LUT查找次数
M_D = 100;           % 总投影方向数
len_skm = 70;        % 水印投影方向数
L = 50;              % 水印比特数
Delta=0.5;


[x, fs] = audioread('A2_0.wav');                             % 读取音频文件,audioread函数会自动完成归一化，所以在计算PSNR峰值信噪比的时候应该用1来计算
fprintf('音频长度: %d 个样本, 采样率: %d Hz\n', length(x), fs);
if size(x, 2) > 1
    x = mean(x, 2);                                          % 对左右声道取平均
end

[ca1, cd1] = dwt(x, 'db1');                                  % 一级DWT变换 - 只处理低频系数
M = length(ca1);  % 低频系数长度

flag = 1;
while ~isempty(flag)
    G = zeros(T, L);                                         % 生成G矩阵，大小为T*L，WLUT=G*m_k，嵌入的水印m_K的长度为L，它是L*1大小的列向量
    for i = 1:T
        j = randi([1, L]);
        G(i, j) = 1;
    end
    sum_G = sum(G);
    flag = find(sum_G < R);                                  %G是每行只有一个1的大小为T*L的矩阵，
end

[S, S_A, A] = S_Gen_Audio(M, M_D, len_skm);                  % 生成投影矩阵
[sk_a, sk_w] = Sk_Gen(G, M_D, A, T, L, len_skm, R);          % 生成密钥
delta = round(Delta * rand(len_skm, 1)) - Delta/2;           % 抖动值δ_j均匀分布在[-Δ/2, Δ/2]

ro = S_A' * ca1;                                             % 生成公共量化项 (服务器端预处理)
ro_q = Delta * round((ro - delta)/Delta) + delta;            % 对投影结果进行量化
em_q = zeros(M, 1);                                          % 计算水印信号
for i = 1:len_skm                                            % 累加每个向量的量化差异
    em_q = em_q + (ro_q(i) - ro(i)) * S_A(:, i);
end

% 加密过程
ELUT = round(sigma_E .* randn(T, 1));                        % 生成加密查找表
c_dwt = ca1;
pad_e = zeros(M_D, 1);;                                      % 实际上可以看成M_D*L大小的矩阵，这个矩阵的每一行都有R个1，即L中有R个1

for j = 1:M_D
    for h = 1:R
        pad_e(j) = pad_e(j) + ELUT(sk_a((j-1)*R + h));
    end
    c_dwt = c_dwt + pad_e(j) * S(:, j);                      % S大小为M*M_D，总投影方向有M_D=100个
end
c_dwt = c_dwt + em_q;

% 联合解密与水印
b_k = randi([0,1],L,1);                                      % 生成L位0-1随机比特
m_k = 2 * b_k - 1;
WLUT = Delta/(4*R) * G * (m_k);                              % 生成个性化水印
DLUT = -ELUT + WLUT;
m_emb_dwt = c_dwt;
pad_d = zeros(M_D, 1);                                       % 可以看成M_D*L大小的矩阵，这个矩阵的每一行都有R个1，即L中有R个1

for j = 1:M_D
    for h = 1:R
        pad_d(j) = pad_d(j) + DLUT(sk_a((j-1)*R + h));
    end
    m_emb_dwt = m_emb_dwt + pad_d(j) * S(:, j);              %S大小为M*M_D，总投影方向有M_D=100个
end

% 计算音频质量指标
x_watermarked = idwt(m_emb_dwt, cd1, 'db1');                 % 逆DWT恢复音频
mse_value = mean((x - x_watermarked).^2);
psnr_value = 20 * log10(1.0 / sqrt(mse_value));
fprintf('水印音频PSNR: %.2f dB\n', psnr_value);

% 水印检测
[ca1_emb, ~] = dwt(x_watermarked, 'db1');                    % 对含水印音频进行DWT变换
Bm = zeros(len_skm, T);                                      % 生成Bm矩阵
for i = 1:len_skm
    for j = 1:R
        Bm(i, sk_w((i-1)*R+j)) = 1;
    end
end

GG = Bm * G;                                                 % Bm*G是列正交矩阵，所以GG*GG'会得出一个乘以R倍的单位矩阵，
wp = S_A' * ca1_emb;
wp_q = (Delta * round((wp - delta)/Delta)) + delta;
e = wp - wp_q;
check=Bm*G*m_k*Delta/(4*R);
% fprintf('提取的水印值θ_k：');fprintf('%.4f ',e);
% fprintf('\n实际的水印值θ_k：');fprintf('%.4f ',check);
% for i=1:L
%     fprintf('%.4f ',e(i)-check(i));
% end
ee = GG' * e;                                                % e=Bm*G*m_k,ee=GG*e,ee=(GG')*GG*m_k=I*m_k,这里的ee就是m_k
b_arb = sign(ee);
b_arb = floor((b_arb + 1)/2);

% 检测结果
if isequal(b_arb, b_k)
    fprintf('无噪声情况下水印提取成功！\n');
else
    fprintf('\n无噪声情况下水印提取失败\n');
end

% 保存水印音频
audiowrite('watermarked_audio.wav', x_watermarked, fs);
fprintf('水印音频已保存为: watermarked.wav\n');
% 播放对比
fprintf("播放原始音频.................\n");
sound(x, fs); pause(length(x)/fs + 2);
fprintf("播放水印音频.................\n");
sound(x_watermarked, fs);