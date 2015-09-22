# The-Great-Gatsby
MATLAB
function res = PermsRep(v,k)
%  PERMSREP Permutations with replacement.
%  
%  PermsRep(v, k) lists all possible ways to permute k elements out of 
%  the vector v, with replacement.

if nargin < 1 || isempty(v)    
    error('v must be non-empty')
else
    n = length(v);
end
 
if nargin < 2 || isempty(k)
    k = n;
end

v = v(:).'; %Ensure v is a row vector
for i = k:-1:1
    tmp = repmat(v, n^(k-i), n^(i-1));
    res(:,i) = tmp(:);
end

function simulatedString = simulateIndep(allowedChar, counts, strLength)
%This function generates a string of randomly chosen symbols. The symbols
%are independent and distributed according to their frequencies in array
%"counts".
%
% INPUTS:
%   allowedChar - an array containing the symbols that can be chosen.
%   counts - the frequencies of the corresponding symbols in "allowedChar".
%   length - the number of characters in the output string.
%
% OUTPUT:
%   simulatedString - the output string of randomly chosen symbols.

%Default length is 1000 characters
if nargin < 3 || strLength < 1
    strLength = 1000;
end

%Normalize the symbol frequencies to create the probability distribution
pdf = counts/sum(counts);

%Create the cumulative distribution
cdf = [0, cumsum(pdf)];

%Monte Carlo simulation of 1st order statistics
rr = rand(1, strLength);
[~, rrLetters] = histc(rr, cdf);
simulatedString = allowedChar(rrLetters);

end

function simulatedString = simulateMarkov(allowedChar, counts, N, strLength)
%This function generates a string of randomly chosen symbols. The symbols
%are distributed according to the N-gram frequencies in array "counts".
%
% INPUTS:
%   allowedChar - an array containing the symbols that can be chosen.
%   counts - the frequencies of the corresponding N-gram symbol 
%            combinations in "allowedChar".
%   N - the length of the letter combinations, whose frequencies are
%       computed in the array counts.
%   length - the number of characters in the output string.
%
% OUTPUT:
%   simulatedString - the output string of randomly chosen symbols.

%Default length is 1000 characters
if nargin < 4 || strLength < 1
    strLength = 1000;
end

%Construct Markov matrix
L = length(allowedChar);
MM = reshape(counts, L, length(counts)/L)';
for ii = 1:length(MM)
    rowSum = sum(MM(ii,:));
    
    %Check if the entire row is empty
    if rowSum
        MM(ii,:) = MM(ii,:)/rowSum;
    end
end

%Simulate the Markov chain
transC = [zeros(size(MM,1),1), cumsum(MM,2)];
states = zeros(1,strLength);    %storage of states
possibleStates = find(sum(MM,2));    %find which states can occur
prevStateAll = PermsRep(1:L, N-1);
initalState = possibleStates(randi(length(possibleStates), 1)); %start at random possible state
states(1:N-1) = prevStateAll(initalState,:);
rr = rand(1,strLength);

%Vector C will aid in helping calculate the row of transC when N > 2 (see below)
C = zeros(N-1, 1);
for ii = 0:(N-2)
    C(end - ii) = L^ii;
end

for ii = N:strLength
    %----------------------------------------------------------------------
    % The columns of matrix transC denote the current state, si, and the
    % rows denote the N-1 previous states, so they are organized as:
    % si-N+1| ...| si-2 | si-1  , or:
    %
    %   row  #1: 1 | ... | 1 | 1
    %   row  #2: 1 | ... | 1 | 2
    %                 ...
    %   row   #L: 1 | ... | 1 | L
    %   row #L+1: 1 | ... | 2 | 1
    %   row #L+2: 1 | ... | 2 | 2
    %                 ...
    %
    % So, we need to convert the N-1 previous states into the row index
    % for this matrix. The formula is as follows:
    %
    % [(si-N+1) - 1]*L^(N-2) + [(si-N+2) - 1]*L^(N-3) + ...
    %                         + [(si-2) - 1]*L^(1) + (si-1)*L^(0)     (*)
    %
    % Take for example, 27 states and N = 4 with s(5:7) = [15, 2, 23].
    %   14*27^(2) + 1*27 + 23 = 10256
    %----------------------------------------------------------------------
    
    %First, fetch the previous N-1 states: si-N+1| ...| si-2 | si-1
    prevState = states((ii-N+1):(ii-1)); 
    %Subtract 1 from all previous N-1 except for the most recent, si-1
    prevState(1:(N-2)) = prevState(1:(N-2)) - 1;
    
    %prevState*C applies formula (*) above to access the correct row of transC
    try
    [~,states(ii)] = histc(rr(ii), transC(prevState*C, :));
    catch
        disp('error')
    end
end
simulatedString = allowedChar(states);

end

