%% Problem 1. 
% in the repository you will find the meannonan.m function we discussed in
% class which produced the mean of a vector of numbers that ignores values
% of NaN or Inf. 

% Part 1. Run the following code:

xx = rand(5); % random 5x5 matrix
xx(3,2) = NaN; %put a NaN in
yy = mean(xx); 
zz = meannonan(xx);

%compare the size of yy and zz. Notice that yy produces a row vector (the
%average of the rows) and has a NaN in the column that contained the NaN whereas 
%zz is a single number (the average of all non-NaN entries in xx). Explain
%this behavior. 

%%  Part 1 Answer
    size(yy)

ans =

     1     5
     
     size(zz)

ans =

     1     1

% The yy variable is a row vector with 5 columns. This variable was
% created by taking the mean of the xx matrix, which finds the mean of the
% all row values within a column, producing a 1 x 5 row vector. Because the mean 
% function can't compute the mean of a set of numbers with a NaN entry, this row vector
% has a NaN entry in the 2nd column position. When the xx variable is
% passed as an argument into the meannonan function, the isnan() function
% is applied on the xx matrix, which returns a logical matrix output
% (suppressed) with a number 1 (indicating TRUE, or that that entry is NaN)
% at the 3rd row, 2nd column position as expected. The meannonan function
% then takes that notin matrix and uses it as an index for the original xx
% matrix (still contains NaN entry), setting it equal to empty brackets [];
% this operation disposes the NaN value and converts the rest of the matrix 
% into a row vector, as the dimensions of the matrix have been manipulated
% and can no longer create a 5 x 5 matrix. The meannonan function finally
% takes the mean of the modified xx row vector, adding up all the columns
% to yield a 1 x 1 dimension matrix, which is stored in the zz variable.
% Hence, the zz variable, which stores the output of the meannonan
% function, has a reasonably different dimension than the yy variable,
% which stores the output of the mean function. 

notin=isnan(xx);
notin

notin =

  5×5 logical array

   0   0   0   0   0
   0   0   0   0   0
   0   1   0   0   0
   0   0   0   0   0
   0   0   0   0   0

xx(notin)=[];

% Part 2. Modify the meannonan code so that it behaves as the mean function
% and produces a row vector where each entry is the average of each column
% and in the column with a NaN, this NaN is ignored. 


%See function code below and in separate file:

function mm=meannonan_rowvector(x)
for ii = 1:length(x)
    for y = 1:size(x,1)
        if isnan(x(y, ii)) == 1
        x(y, ii) = 0;
        end
    end
end
denom = zeros(1,length(x));
for ii = 1:length(x)
    dm = size(x,1);
    for y = 1:size(x,1)
        if x(y,ii) == 0
        dm = dm - 1;
        end
    end 
    denom(ii) = dm;
end
meannonan = zeros(1,5);
for ii = 1:length(x)
    meannonan(ii) = sum(x(:,ii))/(denom(ii));
end
mm=meannonan;

% Demonstration of the utility of the function meannonan_rowvector:

xx = rand(5);
xx(3,2) = NaN
xx =

    0.4942    0.3342    0.5000    0.8594    0.8865
    0.7791    0.6987    0.4799    0.8055    0.0287
    0.7150       NaN    0.9047    0.5767    0.4899
    0.9037    0.0305    0.6099    0.1829    0.1679
    0.8909    0.7441    0.6177    0.2399    0.9787
meannonan_rowvector(xx)
ans =

    0.7566    0.4519    0.6224    0.5329    0.5103

%% Problem 2. ORFs using functions
% In this problem we will use functions to simplify and extend our code from HW1, prob 2 

% Part 1. Fill in the function randdnaseq.m in this repository so that it returns a random sequence
% of length N. 

% See function code below and in separate file: 

function randomSeq = randdnaseq(N)
bases = ['A', 'T', 'G', 'C'];
randseq = bases(randi(4,1,N));
randomSeq = randseq;

% Example of function application: 

randdnaseq(40)

ans =

    'ATTAGATGTTGTCCGTGACCCTGATTAATAGTGGGAATGG'
        
% Part 2. Fill in the function findORF.m in this repository so that takes any dna
% sequence as an input and returns the length of the longest open
% reading frame and the positions of the start and stop codons. 
% Decide what your code should do when no ORF is found and
% implement this. Your function should also work whether the entered dna
% sequence is uppercase, lowercase, or some mixture. The builtin MATLAB functions
% lower and upper could be useful for this. 

% See function findORF() in separate file.

% Part 3. Write another function called probabilityORF that utilizes the functions from 
% Parts 1 and 2. It should take two inputs - a sequence length (N) and an length  of an ORF (N_ORF) and
% returns the probability that that a sequence of length N contains an ORF
% of at least length N_ORF

% See function probabilityORF() in separate file.

% Part4. Write  a final function called plotProbabilityORF.m which takes
% N_ORF as an argument and makes a plot of the probabily of having an
% ORF at least this long as a function of the dnasequence length. Decide how the
% code should determine the lengths of dna sequence to test and implement
% your decision. 

% See function plotProbabilityORF() in separate file.

% Part 5. Write code that uses your function from part 4 to make a single
% plot with separate curves for ORF lengths 100,200,300,and 400. Make sure
% your plot has appropriate axis labels and legend. 

ORF_lengths = [100, 200, 300, 400];
figure;
plotProbabilityORF(ORF_lengths(1)); hold on;
plotProbabilityORF(ORF_lengths(2)); hold on;
plotProbabilityORF(ORF_lengths(3)); hold on;
plotProbabilityORF(ORF_lengths(4)); hold on;
legend('100', '200', '300', '400');
%% Problem 3. Codon translation and optimization

% DNA sequence gets translated into protein through a code known as the
% genetic code. Every sequence of 3 base pairs (a codon) is translated into
% 1 amino acid. The first two columns of the file codons.csv file in this repository 
% give the correspondence between codons and amino acids. 

% Part 1. Fill in the function dna2protein.m so that it takes any
% dnasequence and translates it to protein. The second argument -  frame -
% should take on values of 1,2, or 3 and should refer to whether the
% translation should start from the 1st, 2nd or 3rd base pair (that is,
% which reading frame to use). Make your code returns an error and
% appropriate message if frame isn't 1,2, or 3. 

% See function dna2protein() in separate file.

% Part 2. Write code to turn your protein sequence back into DNA sequence.
% Call your function protein2dna.m
% Notice that there isn't a unique way to do this. For example, there are 4
% different codons that correspond to the amino acid Gly. For a first pass,
% choose one of these codons at random.

% See function protein2dna() in separate file.

% Part 3. The third column of the codons.csv file contains the frequency of
% this codon's use in the human proteome in units of number of appearances per
% thousand codons. Some codons are used more than others. For example,
% for the amino acid Gln, the codon CAG is used nearly 3 times as often as
% CAA. When researchers add DNA to human cells, if it contains CAG rather than
% CAA for Gln, it will be translated  more efficiently. The process of
% taking a protein seqeunce and finding the DNA sequence that will be
% translated most efficiently in a particular organism is called codon
% optimization. Copy your function protein2dna.m from part 2 to a new
% function called protein2dnaOptimized.m that produces a codon-optimized DNA sequence using the 
% information in the third column of codons.csv. 
% In other words, for any amino acid, it always uses the codon that appears
% most frequently in the human proteome. 

% See function protein2dnaOptimized() in separate file. 
