function msg_coded = encode(data_in,G_generator,frozen_positions,free_positions)
[N1,N2]=size(G_generator);
G_frozen=zeros(length(frozen_positions),N1);
G_free=zeros(length(free_positions),N2);
for i=1:1:length(frozen_positions)
    G_frozen(i,:)=G_generator(frozen_positions(i),:);
end
for i=1:1:length(free_positions)
    G_free(i,:)=G_generator(free_positions(i),:);
end
msg_coded=mod(data_in'* G_free,2);
msg_coded=msg_coded';

    