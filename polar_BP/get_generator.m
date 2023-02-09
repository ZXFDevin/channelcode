function G_generator= get_generator(N)
G_generator=zeros(N,N);
n=log2(N);
for row_index=0:1:(N-1)
      b = dec2bin(row_index,n);
      for column_index=0:1:(N-1)
          b_dash = dec2bin(column_index,n);
          value=1;
          for i = 1:1:n 
              value= value*mod((1+b_dash(i)+b(n-i+1)*b_dash(i)),2);
          end
			  G_generator(row_index+1, column_index+1) = value;
      end
end