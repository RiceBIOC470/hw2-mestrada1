function randomSeq = randdnaseq(N)
bases = ['A', 'T', 'G', 'C'];
randseq = bases(randi(4,1,N));
randomSeq = randseq;
% returns a random dna sequence of length N
