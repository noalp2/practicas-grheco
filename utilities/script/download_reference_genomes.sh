wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.basic.annotation.gtf.gz
mv GRCh38.primary_assembly.genome.fa.gz ../../utilities/GRCh38.primary_asembly_reference_genome.fa.gz
mv gencode.v49.transcripts.fa.gz ../../utilities/GCRh38.primary_asembly_reference_transcripts.fa.gz
mv gencode.v49.primary_assembly.basic.annotation.gtf.gz ../../utilities/GCRh38_primary_assembly_reference_annotation.gtf.gz