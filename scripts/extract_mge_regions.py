"""Extract integrative element coordinates from BLAST and tRNA tables."""

"""Extract integrative element coordinates from BLAST and tRNA tables."""

import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from logger import Logger

logger = Logger(name="extract_mge_regions", level=Logger.Level.DEBUG).get_logger()


def parse_blast(blast_file):
    """Parse the BLAST file containing HMM search results.

    The BLAST file (mge_blast.tsv) is expected to have the following columns:
    integrase_id, contig, hit_start, hit_end, pident, length, evalue, bitscore.
    If the file is empty, an empty DataFrame with the specified columns is returned.

    Returns:
        pd.DataFrame: DataFrame containing the best hit per integrase_id.
    """
    if os.stat(blast_file).st_size == 0:
        logger.warning(f"BLAST file {blast_file} is empty. Creating an empty DataFrame.")
        return pd.DataFrame(columns=["integrase_id", "contig", "hit_start", "hit_end", "pident", "length", "evalue", "bitscore"])
    
    df = pd.read_csv(blast_file, sep='\t', header=0, dtype={"contig": str})
    if df.empty:
        logger.warning(f"BLAST file {blast_file} contains no data.")
    else:
        logger.debug("Parsed BLAST DataFrame before filtering:")
        logger.debug(df.head().to_string())

    df_best = df.sort_values("bitscore", ascending=False).drop_duplicates("integrase_id", keep="first")

    logger.debug("BLAST DataFrame after keeping best hits per integrase_id:")
    logger.debug(df_best.head().to_string())

    return df_best


def parse_trna(trna_file):
    """Parse the tRNA file containing integrase-tRNA associations.

    The tRNA file (integrase_trna.tsv) is expected to have the following columns:
    integrase_id, model, integrase_start, integrase_end, integrase_strand,
    contig, contig_length, trna_start, trna_end, trna_strand, tRNA_type, distance.
    If the file is empty, an empty DataFrame with the specified columns is returned.

    Args:
        trna_file (str): Path to the tRNA TSV file.

    Returns:
        pd.DataFrame: DataFrame containing the tRNA data.
    """
    if os.stat(trna_file).st_size == 0:
        logger.warning(f"tRNA file {trna_file} is empty. Creating an empty DataFrame.")
        return pd.DataFrame(columns=[
            "integrase_id", "model", "integrase_start", "integrase_end", "integrase_strand",
            "contig", "contig_length", "trna_start", "trna_end", "trna_strand", "tRNA_type", "distance"
        ])
    
    df = pd.read_csv(trna_file, sep='\t', header=0, dtype={"contig": str})
    if df.empty:
        logger.warning(f"tRNA file {trna_file} contains no data.")
    else:
        logger.debug("Parsed tRNA DataFrame:")
        logger.debug(df.head().to_string())
    return df


def filter_closest_trna(trna_df):
    """Filter the closest tRNA to each integrase based on absolute distance.

    This function finds, for each integrase_id, the tRNA with the smallest absolute
    distance value. If multiple tRNAs exist for the same integrase, only the closest
    one is retained.

    Args:
        trna_df (pd.DataFrame): DataFrame containing tRNA data with a 'distance' column.

    Returns:
        pd.DataFrame: Filtered DataFrame with only one closest tRNA per integrase_id.
    """
    if trna_df.empty:
        return trna_df
    trna_df = trna_df.copy()
    trna_df['abs_distance'] = trna_df['distance'].abs()
    trna_df = trna_df.sort_values('abs_distance').drop_duplicates('integrase_id', keep='first')
    trna_df = trna_df.drop(columns='abs_distance')
    return trna_df


def build_region_data(blast_df, trna_df):
    """Build a DataFrame with region coordinates for MGE extraction.

    For each integrase_id present in the tRNA DataFrame, this function looks up the corresponding
    record in the BLAST DataFrame. The region is calculated based on the tRNA strand:
      - If trna_strand is '+': mge_start = trna_start and mge_end = hit_end.
      - If trna_strand is '-': mge_start = hit_end and mge_end = trna_end.
    If there is no BLAST data for a given integrase_id, that entry is skipped.

    Args:
        blast_df (pd.DataFrame): DataFrame containing BLAST results.
        trna_df (pd.DataFrame): DataFrame containing tRNA data.

    Returns:
        pd.DataFrame: A new DataFrame with columns: integrase_id, contig_id, mge_start, mge_end.
    """
    if blast_df.empty or trna_df.empty:
        logger.warning("One of the input DataFrames is empty. Cannot build regions.")
        return pd.DataFrame(columns=["integrase_id", "contig_id", "mge_start", "mge_end"])
    
    # Build a dictionary for fast lookup of BLAST records by integrase_id.
    blast_dict = {row['integrase_id']: row for _, row in blast_df.iterrows()}
    regions = []
    for _, t_row in trna_df.iterrows():
        integrase_id = t_row['integrase_id']
        if integrase_id not in blast_dict:
            logger.warning(f"No BLAST data for {integrase_id}")
            continue
        blast_row = blast_dict[integrase_id]
        contig_id = blast_row['contig']
        try:
            hit_end = int(blast_row['hit_end'])
            trna_strand = str(t_row['trna_strand']).strip()
            if trna_strand == '+':
                mge_start = int(t_row['trna_start'])
                mge_end = hit_end
            elif trna_strand == '-':
                mge_start = hit_end
                mge_end = int(t_row['trna_end'])
            else:
                logger.error(f"Undefined orientation for {integrase_id}: {trna_strand}")
                continue
            regions.append({
                'integrase_id': integrase_id,
                'contig_id': contig_id,
                'mge_start': mge_start,
                'mge_end': mge_end
            })
        except Exception as e:
            logger.error(f"Error processing {integrase_id}: {e}")
    region_df = pd.DataFrame(regions)
    logger.debug("New DataFrame with region coordinates:")
    logger.debug(region_df.to_string(index=False))
    return region_df


def extract_regions(genome_fasta, region_df, out_fa):
    """Extract regions from the genome based on the coordinates in region_df.

    For each row in region_df, the function extracts the sequence from the genome FASTA file
    corresponding to the region defined by mge_start and mge_end (1-indexed). The FASTA header
    is constructed as: integrase_id:contig_id:mge_start-mge_end.

    Args:
        genome_fasta (str): Path to the genome FASTA file.
        region_df (pd.DataFrame): DataFrame with columns: integrase_id, contig_id, mge_start, mge_end.
        out_fa (str): Path to the output FASTA file.

    Returns:
        None
    """
    records = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    seq_records = []
    for _, row in region_df.iterrows():
        contig_id = str(row['contig_id'])
        mge_start = row['mge_start']
        mge_end = row['mge_end']
        if contig_id in records:
            rec = records[contig_id]
        else:
            rec = None
            for key in records:
                if contig_id in key:
                    rec = records[key]
                    logger.info(f"Using {key} instead of {contig_id}")
                    break
            if rec is None:
                logger.error(f"Contig {contig_id} not found in genome")
                continue
        # Coordinates are 1-indexed.
        seq = rec.seq[mge_start - 1 : mge_end]
        header = f"{row['integrase_id']}:{contig_id}:{mge_start}-{mge_end}"
        seq_records.append(SeqRecord(seq, id=header, description=""))
    with open(out_fa, "w") as f:
        SeqIO.write(seq_records, f, "fasta")
    logger.info(f"Extracted {len(seq_records)} regions. Results saved in {out_fa}")


def main():
    parser = argparse.ArgumentParser(description="Extract MGE regions without merging DataFrames")
    parser.add_argument("--fna", required=True, help="Path to the genome FASTA file")
    parser.add_argument("--blast", required=True, help="Path to the BLAST file (mge_blast.tsv)")
    parser.add_argument("--trna", required=True, help="Path to the tRNA file (integrase_trna.tsv)")
    parser.add_argument("--out_fa", required=True, help="Path to the output FASTA file")
    args = parser.parse_args()
    
    blast_df = parse_blast(args.blast)
    trna_df = parse_trna(args.trna)
    trna_df = filter_closest_trna(trna_df)
    region_df = build_region_data(blast_df, trna_df)
    extract_regions(args.fna, region_df, args.out_fa)


if __name__ == "__main__":
    main()