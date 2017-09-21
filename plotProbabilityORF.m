function output = plotProbabilityORF(N_ORF)
DNAseqlengths = [50, 100, 250, 500, 1000, 2500, 5000, 10000];
probability_vector = zeros(1,8);
for ii=1:8
    probability_vector(ii) = probabilityORF(DNAseqlengths(ii), N_ORF);
end
figure;
plot1=plot(DNAseqlengths, probability_vector);
xlabel('Sequence Length');
ylabel('Probability'); 
output = plot1;

