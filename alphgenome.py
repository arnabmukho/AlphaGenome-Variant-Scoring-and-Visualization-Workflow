from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.interpretation import ism
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
import pandas as pd
from google.colab import userdata
from google.colab import files
import random
from pyfaidx import Fasta

#Create a simulated SNP and an Alphagenome input-ready file
#Download hg38.fa.gz (https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)

bed_file = "test.bed" 
fasta_file = "hg38.fa"      # Indexed FASTA for the correct reference build
out_vcf = "alphagenome_input.vcf"

# ---- LOAD DATA ----
bed = pd.read_csv(bed_file, sep="\t", header=None, names=["CHROM", "START", "END", "GENE"])
fasta = Fasta(fasta_file)
bases = ["A", "C", "G", "T"]

# ---- VCF HEADER ----
vcf_lines = [
    "##fileformat=VCFv4.2",
    "##source=bed_to_vcf_col2.py",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
]

for i, row in bed.iterrows():
    chrom = str(row["CHROM"])
    start0 = int(row["START"])           # BED start (0-based)
    pos1 = start0 + 1                    # VCF is 1-based
    gene = str(row["GENE"]) if not pd.isna(row["GENE"]) else f"gene_{i+1}"

    # Reference base at col2/start
    try:
        ref_base = fasta[chrom][start0].seq.upper()
    except Exception as e:
        print(f"Warning: {chrom}:{pos1} not in FASTA: {e}, skipping.")
        continue

    if ref_base not in bases:
        print(f"Warning: {chrom}:{pos1} non-ACGT ref ({ref_base}), skipping.")
        continue

    alt_base = random.choice([b for b in bases if b != ref_base])

    variant_id = f"{gene}_{chrom}_{pos1}_{ref_base}_{alt_base}"

    vcf_line = f"{chrom}\t{pos1}\t{variant_id}\t{ref_base}\t{alt_base}\t.\t.\t."
    vcf_lines.append(vcf_line)

with open(out_vcf, "w") as f:
    for line in vcf_lines:
        f.write(line + "\n")

print(f"Wrote VCF to {out_vcf}")

# --- Set up API key ---
api_key = os.environ.get("ALPHAGENOME_API_KEY")
if api_key is None:
    raise RuntimeError(
        "Please set the ALPHAGENOME_API_KEY environment variable before running this script.\n"
        "You can do this in your shell with:\nexport ALPHAGENOME_API_KEY=your_real_api_key"
    )

# Initialize the AlphaGenome DNA model using your API key
dna_model = dna_client.create(api_key)

# Uncomment to list output types if you want
print([output.name for output in dna_client.OutputType])

input_file = "alphagenome_input.vcf"
output_file = "alphagenome_variant_scores.tsv"

# ---- Load VCF variants ----
vcf_variants = []
with open(input_file) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 5:
            continue
        chrom, pos, var_id, ref, alt = fields[:5]
        vcf_variants.append({
            'variant_id': var_id,
            'CHROM': chrom,
            'POS': int(pos),
            'REF': ref,
            'ALT': alt
        })
vcf = pd.DataFrame(vcf_variants)
required_cols = ["variant_id", "CHROM", "POS", "REF", "ALT"]
assert all(col in vcf.columns for col in required_cols), f"VCF must contain columns: {required_cols}"

organism = 'human'

sequence_length = '2KB'
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
    f'SEQUENCE_LENGTH_{sequence_length}'
]

score_rna_seq = True
score_cage = True
score_procap = True
score_atac = True
score_dnase = True
score_chip_histone = True
score_chip_tf = True
score_polyadenylation = True
score_splice_sites = True
score_splice_site_usage = True
score_splice_junctions = True

download_predictions = False

organism_map = {
    'human': dna_client.Organism.HOMO_SAPIENS,
    'mouse': dna_client.Organism.MUS_MUSCULUS,
}
organism = organism_map[organism]

scorer_selections = {
    'rna_seq': score_rna_seq,
    'cage': score_cage,
    'procap': score_procap,
    'atac': score_atac,
    'dnase': score_dnase,
    'chip_histone': score_chip_histone,
    'chip_tf': score_chip_tf,
    'polyadenylation': score_polyadenylation,
    'splice_sites': score_splice_sites,
    'splice_site_usage': score_splice_site_usage,
    'splice_junctions': score_splice_junctions,
}

all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
selected_scorers = [
    all_scorers[key]
    for key in all_scorers
    if scorer_selections.get(key.lower(), False)
]

unsupported_scorers = [
    scorer
    for scorer in selected_scorers
    if (
        organism.value
        not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
    )
    | (
        (scorer.requested_output == dna_client.OutputType.PROCAP)
        & (organism == dna_client.Organism.MUS_MUSCULUS)
    )
]
if len(unsupported_scorers) > 0:
    print(
        f'Excluding {unsupported_scorers} scorers as they are not supported for'
        f' {organism}.'
    )
    for unsupported_scorer in unsupported_scorers:
        selected_scorers.remove(unsupported_scorer)

results = []

for _, vcf_row in vcf.iterrows():
    variant = genome.Variant(
        chromosome=str(vcf_row.CHROM),
        position=int(vcf_row.POS),
        reference_bases=vcf_row.REF,
        alternate_bases=vcf_row.ALT,
        name=vcf_row.variant_id,
    )
    interval = variant.reference_interval.resize(sequence_length)

    variant_scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=selected_scorers,
        organism=organism,
    )
    results.append(variant_scores)

df_scores = variant_scorers.tidy_scores(results)

if download_predictions:
    df_scores.to_csv('variant_scores.csv', index=False)
    # files.download('variant_scores.csv')  # Uncomment if using in Colab
else:
    df_scores.to_csv(output_file, sep="\t", index=False)

print(df_scores)

columns = [c for c in df_scores.columns if c != 'ontology_curie']
print(df_scores[(df_scores['ontology_curie'] == 'EFO:0001203')][columns])

gtf = pd.read_feather(
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)

gtf_transcripts = gene_annotation.filter_protein_coding(gtf)
gtf_transcripts = gene_annotation.filter_to_longest_transcript(gtf_transcripts)
transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcripts)

variant = genome.Variant(
    chromosome='chr1',
    position=154220456,
    reference_bases='G',
    alternate_bases='A',
)
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_2KB)

# -- Predict variant effect --
variant_output = dna_model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['EFO:0001203'],  # Example: MCF-7, change as needed
)

# -- Extract relevant transcripts --
longest_transcripts = transcript_extractor.extract(interval)

# -- Plot only the valid intersection of intervals --
plot_interval = variant_output.reference.rna_seq.interval.intersect(interval)

plot_components.plot(
    [
        plot_components.TranscriptAnnotation(longest_transcripts),
        plot_components.OverlaidTracks(
            tdata={
                'REF': variant_output.reference.rna_seq,
                'ALT': variant_output.alternate.rna_seq,
            },
            colors={'REF': 'dimgrey', 'ALT': 'red'},
        ),
    ],
    interval=plot_interval,
    annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
)
plt.show()
