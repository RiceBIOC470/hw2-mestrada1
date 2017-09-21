function proteinseq = dna2protein(dnaseq,frame)
if frame == 1
    a = length(dnaseq);
    num_codons = floor(a/3);
    aminoacidseq = cell([1, num_codons]);
    aminoacidlength = length(aminoacidseq);
    ii = 1;
    for xx = 1:aminoacidlength
        codon = dnaseq(ii:ii+2);
        if strcmp(codon, 'GGG')==1 || strcmp(codon, 'GGA')==1 || strcmp(codon,'GGT')==1 || strcmp(codon,'GGC')==1
                codon = 'Gly';
        elseif strcmp(codon,'GAG')==1 || strcmp(codon,'GAA')==1
                codon = 'Glu';
        elseif strcmp(codon,'GAT')==1 || strcmp(codon,'GAC')==1
                codon = 'Asp';
        elseif strcmp(codon,'GTG')==1 || strcmp(codon,'GTA')==1 || strcmp(codon,'GTT')==1 || strcmp(codon,'GTC')==1
                codon = 'Val';
        elseif strcmp(codon,'GCG')==1 || strcmp(codon,'GCA')==1 || strcmp(codon,'GCT')==1 || strcmp(codon,'GCC')==1
                codon = 'Ala';
        elseif strcmp(codon,'AGG')==1 || strcmp(codon,'AGA')==1 || strcmp(codon,'CGG')==1 || strcmp(codon,'CGA')==1 || strcmp(codon,'CGT')==1 || strcmp(codon,'CGC')==1
                codon = 'Arg';
        elseif strcmp(codon,'AGT')==1 || strcmp(codon,'AGC')==1 || strcmp(codon,'TCG')==1 || strcmp(codon,'TCA')==1 || strcmp(codon,'TCT')==1 || strcmp(codon,'TCC')==1
                codon = 'Ser';
        elseif strcmp(codon,'AAG')==1 || strcmp(codon,'AAA')==1
                codon = 'Lys';
        elseif strcmp(codon,'AAT')==1 || strcmp(codon,'AAC')==1
                codon = 'Asn';
        elseif strcmp(codon,'ATG')==1 
                codon = 'Met';
        elseif strcmp(codon,'ATA')==1 || strcmp(codon,'ATT')==1 || strcmp(codon,'ATC')==1 
                codon = 'Ile';
        elseif strcmp(codon,'ACG')==1 || strcmp(codon,'ACA')==1 || strcmp(codon,'ACT')==1 || strcmp(codon,'ACC')==1 
                codon = 'Thr'; 
        elseif strcmp(codon,'TGG')==1
                codon = 'Trp'; 
        elseif strcmp(codon,'TGT')==1 || strcmp(codon,'TGC')==1 
                codon = 'Cys'; 
        elseif strcmp(codon,'TAT')==1 || strcmp(codon,'TAC')==1
                codon = 'Tyr';
        elseif strcmp(codon,'TTG')==1 || strcmp(codon,'TTA')==1 || strcmp(codon,'CTG')==1 || strcmp(codon,'CTA')==1 || strcmp(codon,'CTT')==1 || strcmp(codon,'CTC')==1
                codon = 'Leu';
        elseif strcmp(codon,'TTT')==1 || strcmp(codon,'TTC')==1
                codon = 'Phe';
        elseif strcmp(codon,'CAG')==1 || strcmp(codon,'CAA')==1 
                codon = 'Gln';
        elseif strcmp(codon,'CAT')==1 || strcmp(codon,'CAC')==1 
                codon = 'His';
        elseif strcmp(codon,'CCG')==1 || strcmp(codon,'CCA')==1 || strcmp(codon,'CCT')==1 || strcmp(codon,'CCC')==1
                codon = 'Pro';
        elseif strcmp(codon,'TGA')==1 || strcmp(codon,'TAG')==1 || strcmp(codon,'TAA')==1
                codon = 'STOP';
        end
        ii = ii + 3;
        aminoacidseq(xx) = {codon};
        if ischar(aminoacidseq{aminoacidlength}) == 1
            break
        end 
    end
elseif frame == 2
    a = length(dnaseq);
    if mod(a,3)==0
    num_codons = floor(a/3)-1;
    elseif mod(a,3)==1
    num_codons = floor(a/3);
    elseif mod(a,3)==2
    num_codons = floor(a/3);
    end
    aminoacidseq = cell([1, num_codons]);
    aminoacidlength = length(aminoacidseq);
    ii = 2;
    for xx = 1:length(aminoacidseq)
        codon = dnaseq(ii:ii+2);
           if strcmp(codon, 'GGG')==1 || strcmp(codon, 'GGA')==1 || strcmp(codon,'GGT')==1 || strcmp(codon,'GGC')==1
                codon = 'Gly';
        elseif strcmp(codon,'GAG')==1 || strcmp(codon,'GAA')==1
                codon = 'Glu';
        elseif strcmp(codon,'GAT')==1 || strcmp(codon,'GAC')==1
                codon = 'Asp';
        elseif strcmp(codon,'GTG')==1 || strcmp(codon,'GTA')==1 || strcmp(codon,'GTT')==1 || strcmp(codon,'GTC')==1
                codon = 'Val';
        elseif strcmp(codon,'GCG')==1 || strcmp(codon,'GCA')==1 || strcmp(codon,'GCT')==1 || strcmp(codon,'GCC')==1
                codon = 'Ala';
        elseif strcmp(codon,'AGG')==1 || strcmp(codon,'AGA')==1 || strcmp(codon,'CGG')==1 || strcmp(codon,'CGA')==1 || strcmp(codon,'CGT')==1 || strcmp(codon,'CGC')==1
                codon = 'Arg';
        elseif strcmp(codon,'AGT')==1 || strcmp(codon,'AGC')==1 || strcmp(codon,'TCG')==1 || strcmp(codon,'TCA')==1 || strcmp(codon,'TCT')==1 || strcmp(codon,'TCC')==1
                codon = 'Ser';
        elseif strcmp(codon,'AAG')==1 || strcmp(codon,'AAA')==1
                codon = 'Lys';
        elseif strcmp(codon,'AAT')==1 || strcmp(codon,'AAC')==1
                codon = 'Asn';
        elseif strcmp(codon,'ATG')==1 
                codon = 'Met';
        elseif strcmp(codon,'ATA')==1 || strcmp(codon,'ATT')==1 || strcmp(codon,'ATC')==1 
                codon = 'Ile';
        elseif strcmp(codon,'ACG')==1 || strcmp(codon,'ACA')==1 || strcmp(codon,'ACT')==1 || strcmp(codon,'ACC')==1 
                codon = 'Thr'; 
        elseif strcmp(codon,'TGG')==1
                codon = 'Trp'; 
        elseif strcmp(codon,'TGT')==1 || strcmp(codon,'TGC')==1 
                codon = 'Cys'; 
        elseif strcmp(codon,'TAT')==1 || strcmp(codon,'TAC')==1
                codon = 'Tyr';
        elseif strcmp(codon,'TTG')==1 || strcmp(codon,'TTA')==1 || strcmp(codon,'CTG')==1 || strcmp(codon,'CTA')==1 || strcmp(codon,'CTT')==1 || strcmp(codon,'CTC')==1
                codon = 'Leu';
        elseif strcmp(codon,'TTT')==1 || strcmp(codon,'TTC')==1
                codon = 'Phe';
        elseif strcmp(codon,'CAG')==1 || strcmp(codon,'CAA')==1 
                codon = 'Gln';
        elseif strcmp(codon,'CAT')==1 || strcmp(codon,'CAC')==1 
                codon = 'His';
        elseif strcmp(codon,'CCG')==1 || strcmp(codon,'CCA')==1 || strcmp(codon,'CCT')==1 || strcmp(codon,'CCC')==1
                codon = 'Pro';
        elseif strcmp(codon,'TGA')==1 || strcmp(codon,'TAG')==1 || strcmp(codon,'TAA')==1
                codon = 'STOP';
        end
        ii = ii + 3;
        aminoacidseq(xx) = {codon};
        if ischar(aminoacidseq{aminoacidlength}) == 1
            break
        end 
    end
elseif frame == 3
    a = length(dnaseq);
    if mod(a,3)==0
    num_codons = floor(a/3)-1;
    elseif mod(a,3)==1
    num_codons = floor(a/3)-1;
    elseif mod(a,3)==2
    num_codons = floor(a/3);
    end
    aminoacidseq = cell([1, num_codons]);
    aminoacidlength = length(aminoacidseq);
    ii = 3;
    for xx = 1:length(aminoacidseq)
        codon = dnaseq(ii:ii+2);
           if strcmp(codon, 'GGG')==1 || strcmp(codon, 'GGA')==1 || strcmp(codon,'GGT')==1 || strcmp(codon,'GGC')==1
                codon = 'Gly';
        elseif strcmp(codon,'GAG')==1 || strcmp(codon,'GAA')==1
                codon = 'Glu';
        elseif strcmp(codon,'GAT')==1 || strcmp(codon,'GAC')==1
                codon = 'Asp';
        elseif strcmp(codon,'GTG')==1 || strcmp(codon,'GTA')==1 || strcmp(codon,'GTT')==1 || strcmp(codon,'GTC')==1
                codon = 'Val';
        elseif strcmp(codon,'GCG')==1 || strcmp(codon,'GCA')==1 || strcmp(codon,'GCT')==1 || strcmp(codon,'GCC')==1
                codon = 'Ala';
        elseif strcmp(codon,'AGG')==1 || strcmp(codon,'AGA')==1 || strcmp(codon,'CGG')==1 || strcmp(codon,'CGA')==1 || strcmp(codon,'CGT')==1 || strcmp(codon,'CGC')==1
                codon = 'Arg';
        elseif strcmp(codon,'AGT')==1 || strcmp(codon,'AGC')==1 || strcmp(codon,'TCG')==1 || strcmp(codon,'TCA')==1 || strcmp(codon,'TCT')==1 || strcmp(codon,'TCC')==1
                codon = 'Ser';
        elseif strcmp(codon,'AAG')==1 || strcmp(codon,'AAA')==1
                codon = 'Lys';
        elseif strcmp(codon,'AAT')==1 || strcmp(codon,'AAC')==1
                codon = 'Asn';
        elseif strcmp(codon,'ATG')==1 
                codon = 'Met';
        elseif strcmp(codon,'ATA')==1 || strcmp(codon,'ATT')==1 || strcmp(codon,'ATC')==1 
                codon = 'Ile';
        elseif strcmp(codon,'ACG')==1 || strcmp(codon,'ACA')==1 || strcmp(codon,'ACT')==1 || strcmp(codon,'ACC')==1 
                codon = 'Thr'; 
        elseif strcmp(codon,'TGG')==1
                codon = 'Trp'; 
        elseif strcmp(codon,'TGT')==1 || strcmp(codon,'TGC')==1 
                codon = 'Cys'; 
        elseif strcmp(codon,'TAT')==1 || strcmp(codon,'TAC')==1
                codon = 'Tyr';
        elseif strcmp(codon,'TTG')==1 || strcmp(codon,'TTA')==1 || strcmp(codon,'CTG')==1 || strcmp(codon,'CTA')==1 || strcmp(codon,'CTT')==1 || strcmp(codon,'CTC')==1
                codon = 'Leu';
        elseif strcmp(codon,'TTT')==1 || strcmp(codon,'TTC')==1
                codon = 'Phe';
        elseif strcmp(codon,'CAG')==1 || strcmp(codon,'CAA')==1 
                codon = 'Gln';
        elseif strcmp(codon,'CAT')==1 || strcmp(codon,'CAC')==1 
                codon = 'His';
        elseif strcmp(codon,'CCG')==1 || strcmp(codon,'CCA')==1 || strcmp(codon,'CCT')==1 || strcmp(codon,'CCC')==1
                codon = 'Pro';
        elseif strcmp(codon,'TGA')==1 || strcmp(codon,'TAG')==1 || strcmp(codon,'TAA')==1
                codon = 'STOP';
        end
        ii = ii + 3;
        aminoacidseq(xx) = {codon};
        if ischar(aminoacidseq{aminoacidlength}) == 1
            break
        end 
    end
else disp('Error: Frame is not 1, 2, or 3')
end  
proteinseq = aminoacidseq;
