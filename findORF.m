function yield = findORF(dnaseq)
for ii = 1:length(dnaseq)
    dnaseq(ii)=upper(dnaseq(ii));
end
output = zeros(1,3);
start_codon = strfind(dnaseq, 'ATG'); 
stop_codons = [strfind(dnaseq, 'TAA') strfind(dnaseq, 'TGA') strfind(dnaseq, 'TAG')];
stop_codon1 = zeros(1, length(start_codon));
for ii = 1:length(start_codon)
    orflengths = stop_codons - start_codon(ii);
    max_length = 1e8;
    orf_index = 0;
    for xx = 1:length(orflengths)
        if orflengths(xx) > 0 && mod(orflengths(xx),3)==0 && orflengths(xx) < max_length
            max_length = orflengths(xx);
            orf_index = xx;
        end
    end
    if orf_index > 0
        stop_codon1(ii) = stop_codons(orf_index);
    else
        stop_codon1(ii) = start_codon(ii);
    end
end
ORFsizes = stop_codon1 - start_codon + 3;
[maxORFlength, ind_max] = max(ORFsizes);
if maxORFlength > 0
   output(1) = maxORFlength;
   output(2) = start_codon(ind_max);
   output(3) = stop_codon1(ind_max);
else
    output(1) = 0;
    output(2) = 0;
    output(3) = 0;
end
yield = output;
% will output max ORF length, start codon position, and stop codon position, in that order. 
                
        
      