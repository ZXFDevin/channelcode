% 输入信道LLR（最右端），初始化最左端LLR（根据冻结比特指示赋值）
% VN代表变量节点，存储LLR信息，包括从右往左L_RL和从左往右L_LR
% 每次信息更新是对VN中的值进行更新
% 对于码长为N的极化码而言，共有log2(N)+1层，第0层是存储输入信道LLR的VN，之后每一层包含N/2个Z型结构
% BP算法中信息传递按右——>左，左——>右的顺序进行
% 每层之间的信息传递发生在“Z”型结构内
% 每个节点存储两个值，分别是从右往左的L值和从左往右的R值
% 重复上述过程直至达到最大迭代次数，对最左侧VN中的值进行比特判决，获得译码结果
% 对于同一Z型结构上左右两侧节点而言
%   从右往左更新时：左侧上节点i中值更新公式  L_RL_左侧(i) = f(L_RL_右侧(i), L_LR_左侧(j)+L_RL_右侧(j))
%                  左侧下节点j中值更新公式  L_RL_左侧(j) = f(L_LR_左侧(i), L_RL_右侧(i))+L_RL右侧(j)
%   从左往右更新时：右侧上节点i中值更新公式  L_LR_右侧(i) = f(L_LR_左侧(i), L_LR_左侧(j)+L_RL_右侧(j))
%                  右侧下节点j中值更新公式  L_LR_右侧(j) = f(L_LR_左侧(i), L_RL_右侧(i))+L_LR左侧(j)
function [decode_msg, tot_iter,flag]=BP_decode(msg_in,N,frozen_array,sigma,max_iter,index_rule,upnode_M, downnode_M, G, tot_iter)
    n=log2(N);
    L = zeros(N,n+1);
    R = zeros(N,n+1);
    flag=0;
     alph = 0.9375;
 %   alph = 0.93252;
    
    msg_in_order = zeros(1,N);
    for i = 1:N
        msg_in_order(i) = msg_in(index_rule(i));
    end
  
    L(:,n+1) = (2*msg_in_order/(sigma^2))';
    for i = 1:N
        if frozen_array(i) == 1
            R(i,1) = 10^100;
        end
    end
    
    for iter=1:max_iter
        % 从右往左更新
        for phy = n:-1:1
            for index=1:N/2
%                   upnode_index = upnode_M(index, phy);
%                   downnode_index = downnode_M(index, phy);
%                   L(upnode_index, phy) = f_function(L(upnode_index,phy+1), L(downnode_index,phy+1)+R(downnode_index,phy));
%                   L(downnode_index, phy) = L(downnode_index, phy+1) + f_function(L(upnode_index,phy+1),R(upnode_index, phy));

                 L(upnode_M(index, phy),phy) = alph*sign(L(upnode_M(index, phy),phy+1))*sign(L(downnode_M(index, phy),phy+1)+R(downnode_M(index, phy),phy))*min(abs(L(upnode_M(index, phy),phy+1)),abs(L(downnode_M(index, phy),phy+1)+R(downnode_M(index, phy),phy)));
                 L(downnode_M(index, phy),phy) = alph*sign(R(upnode_M(index, phy), phy))*sign(L(upnode_M(index, phy),phy+1))*min(abs(R(upnode_M(index, phy), phy)),abs(L(upnode_M(index, phy),phy+1)))+L(downnode_M(index, phy),phy+1);
%                  L(upnode_M(index, phy),phy) = 2*atanh(tanh(L(upnode_M(index, phy),phy+1)/2)*tanh((L(downnode_M(index, phy),phy+1)+R(downnode_M(index, phy),phy))/2));
%                  L(downnode_M(index, phy),phy) = 2*atanh(tanh(R(upnode_M(index, phy), phy)/2)*tanh(L(upnode_M(index, phy),phy+1)/2))+L(downnode_M(index, phy),phy+1);
            end
        end
        % 从左往右更新
        for phy = 1:n
            for index=1:N/2
%                   upnode_index = upnode_M(index, phy);
%                   downnode_index = downnode_M(index, phy);
%                   R(upnode_index, phy+1) = f_function(R(upnode_index, phy), R(downnode_index, phy)+L(downnode_index, phy+1));
%                   R(downnode_index, phy+1) = R(downnode_index, phy)+ f_function(R(upnode_index, phy),L(upnode_index, phy+1));

               R(upnode_M(index, phy),phy+1) = alph*sign(R(upnode_M(index, phy),phy))*sign(L(downnode_M(index, phy),phy+1)+R(downnode_M(index, phy),phy))*min(abs(R(upnode_M(index, phy),phy)),abs(L(downnode_M(index, phy),phy+1)+R(downnode_M(index, phy),phy)));                
               R(downnode_M(index, phy),phy+1) = alph*sign(R(upnode_M(index, phy), phy))*sign(L(upnode_M(index, phy),phy+1))*min(abs(R(upnode_M(index, phy), phy)),abs(L(upnode_M(index, phy),phy+1)))+R(downnode_M(index, phy),phy);
%                  R(upnode_M(index, phy),phy+1) = 2*atanh(tanh(R(upnode_M(index, phy),phy)/2)*tanh((L(downnode_M(index, phy),phy+1)+R(downnode_M(index, phy),phy))/2));
%                  R(downnode_M(index, phy),phy+1) = 2*atanh(tanh(R(upnode_M(index, phy), phy)/2)*tanh(L(upnode_M(index, phy),phy+1)/2))+R(downnode_M(index, phy),phy);
            end
        end
        %停止策略：生成矩阵
        result = L(:,1)+R(:,1);
        codeword = ((L(:,n+1)+R(:,n+1))<0)';
        decode_msg = (result<0)';
        verify = zeros(1,N);
        for i = 1:N
            verify(i) = codeword(index_rule(i));
        end
        if mod(decode_msg*G,2) == verify
            flag=1;
            tot_iter = tot_iter+iter;
            break;
       end
    end
end

    
    



