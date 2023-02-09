% �����ŵ�LLR�����Ҷˣ�����ʼ�������LLR�����ݶ������ָʾ��ֵ��
% VN��������ڵ㣬�洢LLR��Ϣ��������������L_RL�ʹ�������L_LR
% ÿ����Ϣ�����Ƕ�VN�е�ֵ���и���
% �����볤ΪN�ļ�������ԣ�����log2(N)+1�㣬��0���Ǵ洢�����ŵ�LLR��VN��֮��ÿһ�����N/2��Z�ͽṹ
% BP�㷨����Ϣ���ݰ��ҡ���>���󡪡�>�ҵ�˳�����
% ÿ��֮�����Ϣ���ݷ����ڡ�Z���ͽṹ��
% ÿ���ڵ�洢����ֵ���ֱ��Ǵ��������Lֵ�ʹ������ҵ�Rֵ
% �ظ���������ֱ���ﵽ�������������������VN�е�ֵ���б����о������������
% ����ͬһZ�ͽṹ����������ڵ����
%   �����������ʱ������Ͻڵ�i��ֵ���¹�ʽ  L_RL_���(i) = f(L_RL_�Ҳ�(i), L_LR_���(j)+L_RL_�Ҳ�(j))
%                  ����½ڵ�j��ֵ���¹�ʽ  L_RL_���(j) = f(L_LR_���(i), L_RL_�Ҳ�(i))+L_RL�Ҳ�(j)
%   �������Ҹ���ʱ���Ҳ��Ͻڵ�i��ֵ���¹�ʽ  L_LR_�Ҳ�(i) = f(L_LR_���(i), L_LR_���(j)+L_RL_�Ҳ�(j))
%                  �Ҳ��½ڵ�j��ֵ���¹�ʽ  L_LR_�Ҳ�(j) = f(L_LR_���(i), L_RL_�Ҳ�(i))+L_LR���(j)
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
        % �����������
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
        % �������Ҹ���
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
        %ֹͣ���ԣ����ɾ���
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

    
    



