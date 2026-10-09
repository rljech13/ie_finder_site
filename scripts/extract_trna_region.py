"""Extract tRNA sequences paired with integrases for downstream BLAST."""

import argparse

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from genome import contig_record, read_genome
from logger import get_logger
from v3ps_filters import closest_trna_rows

logger = get_logger("extract_trna_region")


def extract_trna_sequence(genome_fasta, trna_table_path, out_fa):
    """Write the tRNA closest to each integrase as a BLAST query.

    The sequence is taken on the tRNA strand: reverse-complemented for a
    ``-`` tRNA. The FASTA id is ``integrase_id:contig:start-end:strand``.

    Args:
        genome_fasta (str): Path to the genome FASTA file.
        trna_table_path (str): Path to ``integrase_trna.tsv``.
        out_fa (str): Path to the output FASTA file.
    """
    sequences = read_genome(genome_fasta)
    # Read contig names as text: numeric names like "1" would otherwise become int.
    trna_df = pd.read_csv(trna_table_path, sep="\t", dtype={"contig": str})
    if trna_df.empty:
        logger.info("tRNA table is empty, creating an empty output file.")
        open(out_fa, "w").close()
        return

    records = []
    for _, row in closest_trna_rows(trna_df).iterrows():
        integrase_id = row["integrase_id"]
        contig = str(row["contig"]).strip()
        start = int(row["trna_start"])
        end = int(row["trna_end"])
        strand = str(row["trna_strand"]).strip()

        if start < 1:
            logger.warning(f"Start coordinate {start} is less than 1 for {integrase_id}")
            continue
        rec = contig_record(sequences, contig)
        if end > len(rec.seq):
            logger.warning(f"End coordinate {end} exceeds length of {contig} (length {len(rec.seq)}) for {integrase_id}")
            continue

        seq = rec.seq[start - 1:end]
        if strand == "-":
            seq = seq.reverse_complement()
        if not seq:
            logger.warning(f"Empty sequence for {integrase_id}|{contig}|{start}-{end}|{strand}, skipping")
            continue

        header = f"{integrase_id}:{contig}:{start}-{end}:{strand}"
        records.append(SeqRecord(seq, id=header, description=""))

    with open(out_fa, "w") as f:
        SeqIO.write(records, f, "fasta")
    logger.info(f"Wrote {len(records)} sequences to {out_fa}")


def main():
    """Parse command-line arguments and extract the closest tRNA sequences.

    This function reads the genome FASTA file, tRNA table, and output FASTA file paths from the
    command line, then calls `extract_trna_sequence` to perform the extraction.

    Returns:
        None
    """
    parser = argparse.ArgumentParser(description="Extract tRNA sequences based on given coordinates")
    parser.add_argument("--ffn", required=True, help="Path to the genome FASTA file")
    parser.add_argument("--trnas", required=True, help="Path to the tRNA table (TSV)")
    parser.add_argument("--out_fa", required=True, help="Path to the output FASTA file")
    args = parser.parse_args()
    extract_trna_sequence(args.ffn, args.trnas, args.out_fa)


if __name__ == "__main__":
    main()