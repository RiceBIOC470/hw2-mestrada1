function output = protein2dna(protein_seq)
gly_codons = {'GGG', 'GGA', 'GGT', 'GGC'}; 
glu_codons = {'GAG', 'GAA'};
asp_codons = {'GAT', 'GAC'};
val_codons = {'GTG', 'GTA', 'GTT', 'GTC'};
ala_codons = {'GCG', 'GCA', 'GCT', 'GCC'};
arg_codons = {'AGG', 'AGA','CGG', 'CGA', 'CGT', 'CGC'};
ser_codons = {'AGT', 'AGC', 'TCG', 'TCA', 'TCT', 'TCC'};
lys_codons = {'AAG', 'AAA'};
asn_codons = {'AAT', 'AAC'};
met_codons = {'ATG'};
ile_codons = {'ATA', 'ATT', 'ATC'};
thr_codons = {'ACG', 'ACA', 'ACT', 'ACC'};
trp_codons = {'TGG'};
cys_codons = {'TGT', 'TGC'};
tyr_codons = {'TAT', 'TAC'};
leu_codons = {'TTG', 'TTA', 'CTG', 'CTA', 'CTT', 'CTC'};
phe_codons = {'TTT', 'TTC'};
gln_codons = {'CAG', 'CAA'};
his_codons = {'CAT', 'CAC'};
pro_codons = {'CCG', 'CCA', 'CCT', 'CCC'};
stop_codons = {'TAG', 'TAA', 'TGA'};
dnaseq = cell([1, length(protein_seq)]);
for ii = 1:length(protein_seq)
    if strcmp(protein_seq(ii), 'Gly')==1
        dnaoutput = gly_codons(randi(numel(1:4)));
    elseif strcmp(protein_seq(ii), 'Glu')==1
        dnaoutput = glu_codons(randi(numel(1:2)));
    elseif strcmp(protein_seq(ii), 'Asp')==1
        dnaoutput = asp_codons(randi(numel(1:2)));
    elseif strcmp(protein_seq(ii), 'Val')==1
        dnaoutput = val_codons(randi(numel(1:4)));
    elseif strcmp(protein_seq(ii), 'Ala')==1
        dnaoutput = ala_codons(randi(numel(1:4)));
    elseif strcmp(protein_seq(ii), 'Arg')==1
        dnaoutput = arg_codons(randi(numel(1:6)));
    elseif strcmp(protein_seq(ii), 'Ser')==1
        dnaoutput = ser_codons(randi(numel(1:6)));
    elseif strcmp(protein_seq(ii), 'Lys')==1
        dnaoutput = lys_codons(randi(numel(1:2)));
    elseif strcmp(protein_seq(ii), 'Asn')==1
        dnaoutput = asn_codons(randi(numel(1:2)));
    elseif strcmp(protein_seq(ii), 'Met')==1
        dnaoutput = met_codons(randi(numel(1)));
    elseif strcmp(protein_seq(ii), 'Ile')==1
        dnaoutput = ile_codons(randi(numel(1:3)));
    elseif strcmp(protein_seq(ii), 'Thr')==1
        dnaoutput = thr_codons(randi(numel(1:4)));
    elseif strcmp(protein_seq(ii), 'Trp')==1
        dnaoutput = trp_codons(randi(numel(1)));
    elseif strcmp(protein_seq(ii), 'Cys')==1
        dnaoutput = cys_codons(randi(numel(1:2)));
    elseif strcmp(protein_seq(ii), 'Tyr')==1
        dnaoutput = tyr_codons(randi(numel(1:2)));
    elseif strcmp(protein_seq(ii), 'Leu')==1
        dnaoutput = leu_codons(randi(numel(1:6)));
    elseif strcmp(protein_seq(ii), 'Phe')==1
        dnaoutput = phe_codons(randi(numel(1:2)));
    elseif strcmp(protein_seq(ii), 'Gln')==1
        dnaoutput = gln_codons(randi(numel(1:2)));
    elseif strcmp(protein_seq(ii), 'His')==1
        dnaoutput = his_codons(randi(numel(1:2)));
    elseif strcmp(protein_seq(ii), 'Pro')==1
        dnaoutput = pro_codons(randi(numel(1:4)));
    elseif strcmp(protein_seq(ii), 'STOP')==1
        dnaoutput = stop_codons(randi(numel(1:3)));
    end
   dnaseq(ii) = dnaoutput; 
end
output = dnaseq;