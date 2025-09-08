curl -L https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz -o utilities/reference_genomes_dna/GRCh38.primary_assembly.genome.gtf.gz
gunzip utilities/reference_genomes_dna/GRCh38.primary_assembly.genome.gtf.gz
curl -L https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz -o utilities/reference_genomes_rna/gencode.v44.annotation.gtf.gz
gunzip utilities/reference_genomes_rna/gencode.v44.annotation.gtf.gz
mv utilities/reference_genomes_dna/GRCh38.primary_assembly.genome.gtf utilities/reference_genomes_rna

