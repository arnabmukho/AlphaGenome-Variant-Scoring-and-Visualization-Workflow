# AlphaGenome Workflow: Scoring Simulated SNP Variants from BED Files

This repository provides a workflow for scoring simulated SNP variants using the AlphaGenome API, starting from a BED file of variants. The workflow reads SNP coordinates from a BED file, generates variant objects, computes functional impact scores (across multiple regulatory tracks), and outputs results in a tidy tabular format suitable for downstream analysis. Optionally, it can visualize the regulatory impact of selected variants.

---

## Features

- **Reads SNP variants from a BED file** (chrom, start, end, name, ref, alt)
- **Scores variants using AlphaGenome deep learning models**
- **Supports diverse functional genomics tracks** (RNA-seq, ATAC-seq, ChIP, etc.)
- **Outputs a tidy TSV table of variant scores**
- **Visualizes variant impact at the transcript and regulatory track level**
- **Secure API authentication via environment variable**

---

## Requirements

- Python 3.8+
- [AlphaGenome Python SDK](https://pypi.org/project/alphagenome/)
- pandas
- matplotlib
- (for feather-format GTFs) pyarrow

Install dependencies with:
```bash
pip install alphagenome pandas matplotlib pyarrow
```

---

## Input BED File Format

The input BED file should contain one SNP per line, with the following columns:
```
chrom   start   end   name   ref   alt
chr1    12345   12346 snp1   A     G
chr2    67890   67891 snp2   C     T
...
```
- `start` is 0-based, `end` is 1-based (standard BED).
- Only SNPs (single nucleotide polymorphisms) are supported.

---

## Usage

1. **Set your AlphaGenome API key**  
   ```bash
   export ALPHAGENOME_API_KEY="your_real_api_key"
   ```

2. **Prepare your input BED file**  
   Place your BED file (e.g., `simulated_snps.bed`) in the repository directory.

3. **Edit the script if needed**  
   - Set `input_file = "simulated_snps.bed"` in `alphgenome.py`
   - Set `output_file = "variant_scores.tsv"` (or your preferred name)

4. **Run the workflow**
   ```bash
   python3 alphgenome.py
   ```

   The script will:
   - Parse SNPs from the BED file
   - Score each SNP using AlphaGenome
   - Output a TSV table of results to `variant_scores.tsv`

5. **(Optional) Visualize variant effects**
   - The script can plot results for a selected variant and cell type/cell line.

---

## Script Overview

- Reads SNPs from a BED file and converts to AlphaGenome variant objects.
- Authenticates with AlphaGenome using your API key.
- Scores each SNP for a suite of regulatory features.
- Outputs a tidy DataFrame (TSV) with all scores.
- Optionally visualizes predicted impact for selected variants.

---

## Plotting

The visualization section selects a SNP for a specified ontology (e.g., cell type like `EFO:0001203` for MCF-7), and plots:
- Reference and alternate predicted tracks (e.g., RNA-seq)
- Transcript annotation overlapping the SNP

---

## Customization

- **Change input/output file names** by editing `input_file` and `output_file` variables in the script.
- **Select scoring tracks** by toggling the `score_*` booleans.
- **Choose ontology/cell type for plotting** by editing the `ontology` variable in the visualization section.

---

## Notes

- **Protobuf Warning**: You may see a warning about protobuf version mismatches. This does not usually affect results, but you can upgrade `alphagenome` or protobuf for best compatibility.
- **Security**: Never commit your API key. Always use environment variables for credentials.

---

## Example Output

- `variant_scores.tsv` — Table of scores for all simulated SNPs.
- Plots — Visualizations of regulatory impact for selected SNPs.

---

## License

See [LICENSE](LICENSE) for details.

---

## Citation

If you use AlphaGenome in your research, please cite the relevant publications and the AlphaGenome platform.

---

For help or contributions, open an issue or pull request.
