function [S, S_A, A] = S_Gen_Audio(M, M_D, len_skm)
    % 输入：
    %   M_D：标准正交向量的总数量
    %   len_skm：需要选取的子集向量数量
    % 输出：
    %   S：M×M_D的标准正交向量集
    %   S_A：从S中选取的len_skm个标准正交向量构成的矩阵
    %   A：选取的向量在S中的索引

    % 随机生成M×(M_D+1)矩阵，通过orth函数正交化
    S = rand(M, M_D+1);
    S = orth(S);
    S(:, 1) = [];            % 移除第一列，保留M_D个向量
    
    % 随机选取len_skm个向量的索引
    M_a = randperm(M_D);     % M_a中是将从1到M_D打乱的索引，M_D=100
    A = M_a(1:len_skm);      % 选择前面len_skm=70个，由于前面进行了随机打乱，所以相当于随机选择所有100投影方向中的70个作为水印投影方向
    A = sort(A);
    
    % 从S中提取选中索引对应的向量，构成S_A
    S_A = zeros(M, len_skm);
    for i = 1:len_skm
        S_A(:, i) = S(:, A(i));
    end
end