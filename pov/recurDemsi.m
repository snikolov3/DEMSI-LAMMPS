function match = recurDemsi(start, target, candidates)

ids = [];
ids = [ids; candidates(candidates(:,1)==start,:); candidates(candidates(:,2)==start,end:-1:1)];
next = ids(ids ~= start);
if next == target
    match = 1;
    return;
else
    candidates(ids(ids ~= start);
    candidates(start, :) = [];
    match = 0;
    
end

% start = starting id in candidates array
% target = target end id terminating triplet 
% for 1 - 2 bond start would be the index for 1 - 2 bond
% and target would be 2 since the triplet would be 1 - X - 2
% 
% In the example below X = 8
% candidates = bond array from base.* file
%
%     1 - 8 - 2  is a triplet
%     7 - 14 - 3 is not a triplet  
%
%
%           1*                                        2*
%         
%      7^         8*                            8*         3^
%      
%  14^    25   31    9                       31     9   14^   32
  
  
  
%  1      7
%  1      8
%  7      14
%  7      25
%  8      31
%  8      9
%  2      8
%  2      3
%  8      31
%  8      9
%  3      14
%  3      32
  
%  ids1 = [1 7; 1 8]
  
 %                               8
 %                             ^   ^
 %                            ^     ^
 %                           ^       ^
 %                          ^         ^
 %                         ^           ^
 %                        ^             ^
 %                       ^               ^
 %                      ^                 ^
 %                     ^                   ^
 %                    1---------------------2