"""Extract tRNA sequences paired with integrases for downstream BLAST."""

"""Extract tRNA sequences paired with integrases for downstream BLAST."""

import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from logger import Logger

logger = Logger(name="extract_trna_region", level=Logger.Level.INFO).get_logger()


def normalize_id(seq_id):
    """Normalize a sequence ID by stripping whitespace and converting to lowercase.

    Args:
        seq_id (str): The sequence identifier.

    Returns:
        str: The normalized sequence identifier.
    """
    return str(seq_id).strip().lower()


def get_sequence_for_contig(sequences, contig):
    """Return the SeqRecord for a given contig.

    If an exact match is not found, searches for a key containing the contig as a substring.

    Args:
        sequences (dict): Dictionary of SeqRecords keyed by normalized IDs.
        contig (str): The contig identifier to search for.

    Returns:
        SeqRecord or None: The matching SeqRecord, or None if not found.
    """
    norm_contig = normalize_id(contig)
    # Look for an exact match.
    if norm_contig in sequences:
        return sequences[norm_contig]
    # If not found, search by substring.
    for key, rec in sequences.items():
        if norm_contig in key:
            logger.info(f"Using {key} for contig {contig}")
            return rec
    logger.error(f"Sequence not found for contig {contig}. Available keys: {list(sequences.keys())}")
    return None


def extract_trna_sequence(genome_fasta, trna_table_path, out_fa):
    """Extract the tRNA sequence closest to each integrase from the genome.

    This function loads genome sequences from a FASTA file and a tRNA table (TSV) containing
    tRNA coordinates. It groups the table by 'integrase_id' and selects only the row with the smallest
    'distance' for each integrase. For tRNAs on the '+' strand, the sequence is extracted as:
        record.seq[trna_start-1:trna_end]
    For tRNAs on the '-' strand, the extracted sequence is reverse complemented.
    The FASTA header will include the integrase_id, contig, coordinate range, and strand.

    Args:
        genome_fasta (str): Path to the genome FASTA file.
        trna_table_path (str): Path to the tRNA coordinates TSV file.
        out_fa (str): Path to the output FASTA file.

    Returns:
        None
    """
    # Load genome sequences and normalize IDs.
    sequences = {normalize_id(rec.id): rec for rec in SeqIO.parse(genome_fasta, "fasta")}
    # Load the tRNA table.
    # Read contig names as text: numeric names like "1" would otherwise become int.
    trna_df = pd.read_csv(trna_table_path, sep="\t", dtype={"contig": str})
    if trna_df.empty:
        logger.info("tRNA table is empty, creating an empty output file.")
        open(out_fa, "w").close()
        return

    # If a 'distance' column exists, keep only the row with the smallest distance per integrase.
    if "distance" in trna_df.columns:
        trna_df = trna_df.loc[trna_df.groupby("integrase_id")["distance"].idxmin()]

    records = []
    for index, row in trna_df.iterrows():
        try:
            integrase_id = row["integrase_id"]
            contig = str(row["contig"]).strip().rstrip(',')
            start = int(row["trna_start"])
            end = int(row["trna_end"])
            strand = str(row["trna_strand"]).strip()
        except Exception as e:
            logger.error(f"Error processing row {row}: {e}")
            continue

        # Check coordinate validity.
        if start < 1:
            logger.warning(f"Start coordinate {start} is less than 1 for {integrase_id}")
            continue

        # Find the corresponding genome sequence for the contig.
        contig_key = normalize_id(contig)
        if contig_key not in sequences:
            # Search for a substring match.
            found = False
            for key in sequences:
                if contig_key in key:
                    contig_key = key
                    found = True
                    break
            if not found:
                logger.error(f"Contig {contig} not found in genome")
                continue

        rec = sequences[contig_key]
        if end > len(rec.seq):
            logger.warning(f"End coordinate {end} exceeds length of {contig} (length {len(rec.seq)}) for {integrase_id}")
            continue

        # Extract the region.
        seq = rec.seq[start - 1:end]
        if strand == "-":
            seq = seq.reverse_complement()

        if not seq or len(seq) == 0:
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