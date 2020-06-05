function parts=partitions(order,next,parts)                                % ‰»Îorder–≈œ¢°£ Next = 0, parts = {{}}
if nargin < 3
   next=0;
   parts={{}};
end
if next < order                                                            %
   next=next + 1;                                                          % next + 1 
   N=numel(parts);                                                         % N = 0
   m=N;                                                                    % m = N
   for n=1:N                                                               % 
      P=numel(parts{n});
      for p=1:P
         m=m + 1;
         parts{m}=parts{n};
         parts{m}{p}(end+1)=next;
      end
      parts{n}{P+1}=next;
   end
   parts=partitions(order,next,parts);
end
