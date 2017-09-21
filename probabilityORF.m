function output=probabilityORF(N, N_ORF)
longest_ORF = zeros(1,1000);
longest_ORFlogical = zeros(1,1000);
for ii=1:1000
    rand_seq = randdnaseq(N);
    rawORFdata=findORF(rand_seq);
    rawORFdata=rawORFdata(1);
    longest_ORF = rawORFdata;
    if longest_ORF >= N_ORF
       longest_ORFlogical(ii)=1;
    end
end
x = sum(longest_ORFlogical);
output=x/1000;

