"""Run HMMER integrase screening and merge hits with ORF coordinates."""

import os
import re
import subprocess
import argparse
from datetime import datetime
from pathlib import Path
from logger import get_logger

logger = get_logger("hmm_search")


HMMSCAN_LONG_SEQ_MSG = "Target sequence length > 100K"


def append_skip_log(skip_log, sample, reason, details=""):
    """Append one skip event to a tab-separated log file.

    Args:
        skip_log: Path to the skip log file, or None to disable logging.
        sample: Sample identifier associated with the skip event.
        reason: Short machine-readable skip reason code.
        details: Optional free-text details.
    """
    if not skip_log:
        return
    Path(skip_log).parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().isoformat()
    with open(skip_log, "a") as log_f:
        log_f.write(f"{timestamp}\t{sample}\t{reason}\t{details}\n")


def run_hmmscan(faa_path, hmm_path, output_tbl, sample=None, skip_log=None):
    """Run hmmscan against a combined Pfam HMM database.

    Args:
        faa_path: Path to the protein FASTA file (``orfs.faa``).
        hmm_path: Path to the pressed combined HMM file.
        output_tbl: Path for hmmscan tabular output (tblout format).
        sample: Optional sample name used in log messages.
        skip_log: Optional path for recording skipped samples.

    Returns:
        Tuple ``(success, skip_reason)`` where ``success`` is True when hmmscan
        completed normally and ``skip_reason`` is a reason code when skipped.

    Raises:
        subprocess.CalledProcessError: When hmmscan fails for a reason other
            than exceeding the 100 K amino-acid target limit.
    """
    logger.info(f"Executing hmmscan for {faa_path} using {hmm_path}...")
    cmd = [
        "hmmscan",
        "--tblout", output_tbl,
        hmm_path,
        faa_path
    ]
    result = subprocess.run(
        cmd,
        text=True,
        capture_output=True
    )
    if result.returncode != 0:
        stderr = result.stderr or ""
        if HMMSCAN_LONG_SEQ_MSG in stderr:
            reason = "hmmscan_target_gt_100k"
            logger.warning(
                f"hmmscan aborted for {sample or faa_path}: "
                f"{HMMSCAN_LONG_SEQ_MSG}. Skipping sample."
            )
            append_skip_log(skip_log, sample or Path(faa_path).stem, reason, HMMSCAN_LONG_SEQ_MSG)
            # Write note to tblout to avoid downstream FileNotFoundError
            Path(output_tbl).write_text(
                f"# Skipped hmmscan for {sample or faa_path}: {HMMSCAN_LONG_SEQ_MSG}\n"
            )
            return False, reason
        else:
            logger.error(
                f"hmmscan failed for {sample or faa_path} with return code {result.returncode}."
            )
            logger.error(stderr)
            result.check_returncode()
    logger.info(f"Results saved in {output_tbl}")
    return True, None


def parse_tblout(tbl_path):
    """Parse the tblout file and create a mapping from ORF ID to model accession.

    The tblout file (e.g., integrase_hits.txt) is expected to have lines in the format:
      [Some header text]
      Phage_integrase      PF00589.27 JBKBIM010000027.1_8  - 4.5e-35 108.1 ...
    where the query name (ORF ID) is the third column and the model accession is the second column.

    Args:
        tbl_path (str): Path to the tblout file.

    Returns:
        dict: A dictionary mapping ORF IDs (str) to model accessions (str).
    """
    mapping = {}
    with open(tbl_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            orf_id = parts[2]
            model_acc = parts[1]
            mapping[orf_id] = model_acc
    return mapping


def parse_faa(faa_path):
    """Parse the FAA file to extract ORF coordinates.

    The FAA file header is expected to have the format:
      >JBKBIM010000001.1_1 # 3 # 1430 # -1 # ID=1_1;partial=10;...
    The function returns a dictionary mapping the ORF ID to a tuple of (start, end, strand).

    Args:
        faa_path (str): Path to the FAA file.

    Returns:
        dict: A dictionary mapping ORF IDs (str) to tuples (start (int), end (int), strand (str)).
    """
    orf_coords = {}
    with open(faa_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                parts = header.split(' # ')
                if len(parts) < 4:
                    logger.error(f"Incorrect FAA header format: {line.strip()}")
                    continue
                orf_id = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    logger.error(f"Error converting coordinates in FAA: {line.strip()}")
                    continue
                strand = parts[3]
                orf_coords[orf_id] = (start, end, strand)
    return orf_coords


def parse_gff(gff_path):
    """Parse the GFF file to extract contig lengths.

    The GFF file is expected to contain a line with the following format:
      # Sequence Data: seqnum=1;seqlen=128375;seqhdr="JBKBIM010000027.1 MAG: Thermus sp. isolate ..."
    This function extracts the sequence length (seqlen) and the contig ID (the first token in seqhdr).

    Args:
        gff_path (str): Path to the GFF file.

    Returns:
        dict: A dictionary mapping contig IDs (str) to their lengths (int).
    """
    contig_lengths = {}
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("# Sequence Data:"):
                m = re.search(r'seqlen=(\d+);seqhdr="([^"]+)"', line)
                if m:
                    seqlen = int(m.group(1))
                    seqhdr = m.group(2)
                    contig_id = seqhdr.split()[0]
                    contig_lengths[contig_id] = seqlen
                else:
                    logger.error(f"Incorrect GFF line: {line.strip()}")
    return contig_lengths


def write_outputs(mapping, faa_coords, contig_lengths, summary_file, orfs_file):
    """Write the integrase summary and ORF coordinate files.

    The summary file contains the following columns:
      ORF ID, model accession, start, end, strand, contig ID, contig length.
    The ORFs file contains the following columns:
      ORF ID, start, end.

    Args:
        mapping (dict): Mapping from ORF ID to model accession.
        faa_coords (dict): Mapping from ORF ID to (start, end, strand).
        contig_lengths (dict): Mapping from contig ID to contig length.
        summary_file (str): Path to the output summary file.
        orfs_file (str): Path to the output ORFs file.
    """
    with open(summary_file, 'w') as summ_f, open(orfs_file, 'w') as orfs_f:
        # Write header for summary file
        summ_f.write("orf_id\tmodel_accession\tstart\tend\tstrand\tcontig_id\tcontig_length\n")
        # Write header for ORFs file
        orfs_f.write("orf_id\tstart\tend\n")
        for orf_id, model in mapping.items():
            if orf_id not in faa_coords:
                logger.error(f"ORF {orf_id} not found in FAA")
                continue
            start, end, strand = faa_coords[orf_id]
            # Contig ID is the ORF ID with the trailing _<orf_index> suffix removed.
            contig_id = '_'.join(orf_id.split('_')[:-1])
            if contig_id not in contig_lengths:
                logger.error(f"Contig {contig_id} not found in GFF")
                continue
            contig_len = contig_lengths[contig_id]
            summ_f.write(f"{orf_id}\t{model}\t{start}\t{end}\t{strand}\t{contig_id}\t{contig_len}\n")
            orfs_f.write(f"{orf_id}\t{start}\t{end}\n")
    logger.info("Data successfully written to output files.")


def main() -> None:
    """Run hmmscan and write integrase summary tables for one sample."""
    parser = argparse.ArgumentParser(
        description="Match HMM scan results with ORF and contig data"
    )
    parser.add_argument("--faa", required=True, help="Path to FAA file (orfs.faa)")
    parser.add_argument("--gff", required=True, help="Path to GFF file (orfs.gff)")
    parser.add_argument("--out", required=True, help="Path to tblout file (integrase_hits.txt)")
    parser.add_argument("--summary", required=True, help="Output file for integrase summary table")
    parser.add_argument("--orfs", required=True, help="Output file for ORF coordinates")
    parser.add_argument("--combined", required=True, help="Path to combined HMM file")
    parser.add_argument("--sample", required=False, help="Sample name (for logging)")
    parser.add_argument(
        "--skip-log",
        required=False,
        help="File to append skipped samples (tab-separated timestamp, sample, reason, details)"
    )
    args = parser.parse_args()

    # Step 1: Run hmmscan
    success, skip_reason = run_hmmscan(
        args.faa,
        args.combined,
        args.out,
        sample=args.sample,
        skip_log=args.skip_log
    )
    if not success:
        # Produce empty outputs so downstream steps continue
        write_outputs({}, {}, {}, args.summary, args.orfs)
        return
    # Step 2: Parse tblout file to get mapping (ORF ID -> model accession)
    mapping = parse_tblout(args.out)
    # Step 3: Parse FAA file to get ORF coordinates
    faa_coords = parse_faa(args.faa)
    # Step 4: Parse GFF file to get contig lengths
    contig_lengths = parse_gff(args.gff)
    # Step 5: Write output summary and ORF coordinate files
    write_outputs(mapping, faa_coords, contig_lengths, args.summary, args.orfs)


if __name__ == "__main__":
    main()