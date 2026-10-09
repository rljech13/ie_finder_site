"""Identify tRNA genes located near integrase hits on the opposite strand."""

import argparse
import csv
import re
import subprocess
from logger import get_logger

logger = get_logger("trna_proximity")


def parse_integrases(integrase_file):
    """Parse the integrase TSV file into a list of dictionaries.

    Reads a TSV file containing integrase data and converts certain fields to integers.
    It also cleans up the contig_id field. Each row is stored as a dictionary.

    Args:
        integrase_file (str): Path to the integrase TSV file (e.g., integrase_hits_summary.tsv).

    Returns:
        list of dict: A list of dictionaries, each representing an integrase entry.
    """
    logger.info(f"Reading integrases from {integrase_file}")
    integrases = []
    with open(integrase_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                row['start'] = int(row['start'])
                row['end'] = int(row['end'])
                row['contig_length'] = int(row['contig_length'])
                row['contig_id'] = row['contig_id'].strip()
                integrases.append(row)
            except Exception as e:
                logger.error(f"Error processing integrase row {row}: {e}")
    logger.info(f"Parsed {len(integrases)} integrases")
    return integrases


def parse_trna(trna_file):
    """Parse the Aragorn output (plain text) for tRNA entries.

    The text file is expected to have blocks starting with a header line (e.g., 
    '>AE017221.1 Thermus thermophilus HB27, complete genome'), followed by a line
    with the number of genes found (e.g., '51 genes found'). The next lines (equal to
    the number of genes found) should contain tRNA data lines in a format like:
      "1   tRNA-Leu   c[3280,3362]   38   (cag)"
    If coordinates do not contain a prefix (i.e., they start with '['), it is assumed
    that the tRNA is on the forward strand ("+").

    Args:
        trna_file (str): Path to the Aragorn plain text output file.

    Returns:
        list of dict: A list of dictionaries, each containing:
            'contig': (str) contig identifier,
            'start': (int) start coordinate,
            'end': (int) end coordinate,
            'strand': (str) strand ("+" or "-"),
            'tRNA_type': (str) tRNA type.
    """
    trnas = []
    with open(trna_file) as f:
        lines = f.read().splitlines()
    current_contig = None
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('>'):
            # Extract contig ID from header (first token)
            parts = line[1:].split()
            current_contig = parts[0].strip()
            i += 1
            if i >= len(lines):
                break
            # Next line should contain the number of genes found
            gene_line = lines[i].strip()
            m_num = re.match(r'(\d+)', gene_line)
            if m_num:
                num_genes = int(m_num.group(1))
            else:
                num_genes = 0
            i += 1
            # Process num_genes lines
            for _ in range(num_genes):
                if i >= len(lines):
                    break
                gene_line = lines[i].strip()
                fields = gene_line.split()
                if len(fields) >= 3:
                    tRNA_type = fields[1]
                    coord_str = fields[2]
                    # Regular expression: optional prefix (f or c), followed by [start,end]
                    m = re.match(r'(?:([fc]))?\[(\d+),(\d+)\]', coord_str)
                    if m:
                        prefix = m.group(1) if m.group(1) is not None else 'f'
                        start = int(m.group(2))
                        end = int(m.group(3))
                        # If prefix is 'c', assume reverse strand, else forward
                        strand = '-' if prefix == 'c' else '+'
                        trnas.append({
                            'contig': current_contig,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'tRNA_type': tRNA_type
                        })
                    else:
                        logger.error(f"Failed to parse coordinates from: {coord_str}")
                i += 1
        else:
            i += 1
    logger.info(f"Parsed {len(trnas)} tRNAs")
    return trnas


def effective_coord_integrase(integ):
    """Calculate the effective coordinate for an integrase.

    For integrases on the '+' or "1" strand, the effective coordinate is the end.
    For integrases on the '-' or "-1" strand, it is the start.

    Args:
        integ (dict): A dictionary representing an integrase entry, expected to contain
            a 'strand' key and 'start'/'end' coordinates.

    Returns:
        int or None: The effective coordinate, or None if strand is undefined.
    """
    if integ['strand'] in ["1", "+"]:
        return integ['end']
    elif integ['strand'] in ["-1", "-"]:
        return integ['start']
    else:
        return None


def effective_coord_trna(trna):
    """Calculate the effective coordinate for a tRNA.

    For tRNAs on the '-' strand, the effective coordinate is the start.
    For tRNAs on the '+' strand, it is the end.

    Args:
        trna (dict): A dictionary representing a tRNA entry, expected to contain a 'strand'
            key along with 'start' and 'end'.

    Returns:
        int or None: The effective coordinate, or None if strand is undefined.
    """
    if trna['strand'] == '-':
        return trna['start']
    elif trna['strand'] == '+':
        return trna['end']
    else:
        return None


def find_nearby_trnas(integrases, trnas, max_distance=500):
    """Find tRNAs located within a specified distance of integrases.

    Groups tRNA entries by contig and then, for each integrase, finds tRNAs on the
    opposite strand that are within max_distance nucleotides of the integrase's effective
    coordinate.

    Args:
        integrases (list of dict): List of integrase entries.
        trnas (list of dict): List of tRNA entries.
        max_distance (int, optional): Maximum allowed distance (default is 500).

    Returns:
        list of dict: A list of result dictionaries containing integrase and tRNA data,
            including the distance between them.
    """
    logger.info("Searching for tRNAs near integrases...")
    results = []
    # Group tRNA entries by contig
    trna_by_contig = {}
    for t in trnas:
        c = t['contig'].strip()
        trna_by_contig.setdefault(c, []).append(t)
    
    for integ in integrases:
        c = integ['contig_id'].strip()
        if c not in trna_by_contig:
            continue
        eff_integ = effective_coord_integrase(integ)
        if eff_integ is None:
            continue
        for t in trna_by_contig[c]:
            # Only match if integrase and tRNA are on opposite strands
            if integ['strand'] in ["1", "+"] and t['strand'] != '-':
                continue
            if integ['strand'] in ["-1", "-"] and t['strand'] != '+':
                continue
            eff_trna = effective_coord_trna(t)
            if eff_trna is None:
                continue
            if integ['strand'] in ["1", "+"]:
                if eff_trna < eff_integ:
                    continue
                dist = eff_trna - eff_integ
            else:
                if eff_integ < eff_trna:
                    continue
                dist = eff_integ - eff_trna
            if dist <= max_distance:
                results.append({
                    'integrase_id': integ['orf_id'],
                    'model': integ['model_accession'],
                    'integrase_start': integ['start'],
                    'integrase_end': integ['end'],
                    'integrase_strand': integ['strand'],
                    'contig': c,
                    'contig_length': integ['contig_length'],
                    'trna_start': t['start'],
                    'trna_end': t['end'],
                    'trna_strand': t['strand'],
                    'tRNA_type': t.get('tRNA_type', ''),
                    'distance': dist
                })
    logger.info(f"Found {len(results)} tRNAs within {max_distance} nt of integrases")
    return results


def write_results(results, output_file):
    """Write the integrase-tRNA proximity results to a TSV file.

    The output file will have the following columns:
        integrase_id, model, integrase_start, integrase_end, integrase_strand,
        contig, contig_length, trna_start, trna_end, trna_strand, tRNA_type, distance

    Args:
        results (list of dict): List of result dictionaries.
        output_file (str): Path to the output TSV file.
    """
    fields = [
        'integrase_id', 'model', 'integrase_start', 'integrase_end', 'integrase_strand',
        'contig', 'contig_length', 'trna_start', 'trna_end', 'trna_strand', 'tRNA_type', 'distance'
    ]
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for r in results:
            writer.writerow(r)
    logger.info(f"Wrote {len(results)} results to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Find tRNAs on the opposite strand within --max_distance nt of integrases"
    )
    parser.add_argument("--integrases", required=True,
                        help="TSV file with integrase data (e.g., integrase_hits_summary.tsv)")
    parser.add_argument("--trna", required=True,
                        help="Aragorn output file (plain text)")
    parser.add_argument("--output", required=True,
                        help="Output TSV file for integrases and nearby tRNAs")
    parser.add_argument("--max_distance", type=int, default=500,
                        help="Maximum allowed distance (default 500)")
    args = parser.parse_args()

    integrases = parse_integrases(args.integrases)
    trnas = parse_trna(args.trna)
    results = find_nearby_trnas(integrases, trnas, args.max_distance)
    write_results(results, args.output)


if __name__ == "__main__":
    main()