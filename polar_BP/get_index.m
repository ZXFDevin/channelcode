function [upnode_M, downnode_M] = get_index(n)
    upnode_M = zeros(2^(n-1), n);
    downnode_M = zeros(2^(n-1), n);
    for phy = 1:n
        index = [];
        node_length = 2^(phy-1);
        for i = 1:2^(phy-1)
            for j = 0:2^(n-phy)-1
                index = [index, i+2^(phy)*j];
            end
        end
        upnode_M(:,phy) = index;
        downnode_M(:,phy) = index+node_length;
    end
end
