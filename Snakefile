import os
import glob
import sys


configfile: os.path.join(workflow.basedir, "ie_finder_config.yaml")

GENOMES_DIR = config["paths"]["genomes_dir"]
RESULTS_DIR = config["paths"]["results_dir"]
PFAM_LIST = config["pfam_profiles"]
PFAM_CLI_ARGS = " ".join(f"--pfam '{pfam}'" for pfam in PFAM_LIST)
ALL_CANDIDATES = bool(config.get("annotate", {}).get("all_candidates", False))
FINDER = os.path.join(workflow.basedir, "scripts")

sys.path.insert(0, FINDER)
from v3ps_filters import load_thresholds

# Search thresholds live in their own file. Loading it here stops a run with an
# unknown key or a bad value before any step starts.
SEARCH_PARAMS = config["paths"].get(
    "search_params", os.path.join(workflow.basedir, "search_params.yaml")
)
THRESHOLDS = load_thresholds(SEARCH_PARAMS)

COMBINED_HMM = os.path.join(RESULTS_DIR, "combined", "pfam_combined.hmm")


def get_samples():
    fasta_files = glob.glob(os.path.join(GENOMES_DIR, "*.fna"))
    if not fasta_files:
        raise FileNotFoundError(f"No .fna files found in {GENOMES_DIR}")
    return [os.path.splitext(os.path.basename(p))[0] for p in fasta_files]


SAMPLES = get_samples()


def _flat(items):
    """Flatten expand() lists so rule all receives file paths."""
    out = []
    for item in items:
        if isinstance(item, (list, tuple)):
            out.extend(item)
        else:
            out.append(item)
    return out


def all_inputs(wildcards):
    """Genome-coordinate GFF3, GenBank, and the attachment-site table."""
    return _flat([
        COMBINED_HMM,
        expand(os.path.join(RESULTS_DIR, "{sample}", "{sample}.ie.gff3"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "{sample}", "{sample}.ie.gbk"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "{sample}", "attachment_sites_genome.tsv"), sample=SAMPLES),
    ])


rule all:
    input:
        all_inputs,


rule predict_orfs:
    input:
        fna=os.path.join(GENOMES_DIR, "{sample}.fna")
    output:
        gff=os.path.join(RESULTS_DIR, "{sample}", "orfs.gff"),
        ffn=os.path.join(RESULTS_DIR, "{sample}", "orfs.ffn"),
        faa=os.path.join(RESULTS_DIR, "{sample}", "orfs.faa")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "predict_orfs.log")
    params:
        finder=FINDER
    shell:
        """
        python {params.finder}/predict_orfs.py --fna {input.fna} --gff {output.gff} --ffn {output.ffn} --faa {output.faa} > {log} 2>&1
        """

rule build_combined_hmm:
    output:
        COMBINED_HMM
    log:
        os.path.join(RESULTS_DIR, "combined", "build_combined_hmm.log")
    params:
        finder=FINDER,
        pfam_args=PFAM_CLI_ARGS
    shell:
        """
        python {params.finder}/build_combined_hmm.py \
            --output {output} \
            --log {log} \
            {params.pfam_args}
        """

rule hmm_search:
    input:
        faa=os.path.join(RESULTS_DIR, "{sample}", "orfs.faa"),
        gff=os.path.join(RESULTS_DIR, "{sample}", "orfs.gff"),
        hmm=COMBINED_HMM
    output:
        hits=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits.txt"),
        stats=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits_summary.tsv"),
        orfs=os.path.join(RESULTS_DIR, "{sample}", "integrase_orfs.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "hmm_search.log")
    params:
        finder=FINDER,
        skip_log=os.path.join(RESULTS_DIR, "hmmscan_skipped.tsv")
    shell:
        """
        python {params.finder}/hmm_search.py \
            --faa {input.faa} \
            --gff {input.gff} \
            --out {output.hits} \
            --summary {output.stats} \
            --orfs {output.orfs} \
            --combined {input.hmm} \
            --sample {wildcards.sample} \
            --skip-log {params.skip_log} \
            > {log} 2>&1
        """

rule predict_trna:
    input:
        fna=os.path.join(GENOMES_DIR, "{sample}.fna")
    output:
        trna=os.path.join(RESULTS_DIR, "{sample}", "trna.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "predict_trna.log")
    shell:
        """
        aragorn -w -t -o {output.trna} {input.fna} > {log} 2>&1
        """

rule trna_proximity:
    input:
        integrases=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits_summary.tsv"),
        trna=os.path.join(RESULTS_DIR, "{sample}", "trna.tsv")
    output:
        proximity=os.path.join(RESULTS_DIR, "{sample}", "integrase_trna.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "trna_proximity.log")
    params:
        finder=FINDER,
        max_distance=THRESHOLDS.trna_max_distance_bp,
    shell:
        """
        python {params.finder}/annotate_trna_proximity.py --integrases {input.integrases} --trna {input.trna} --output {output.proximity} --max_distance {params.max_distance} > {log} 2>&1
        """

rule extract_trna_region:
    input:
        fasta=os.path.join(GENOMES_DIR, "{sample}.fna"),
        trnas=os.path.join(RESULTS_DIR, "{sample}", "integrase_trna.tsv")
    output:
        out_fa=os.path.join(RESULTS_DIR, "{sample}", "mge_query.fa")
    log:
        os.path.join(RESULTS_DIR, "{sample}", "extract_trna_region.log")
    params:
        finder=FINDER
    shell:
        """
        python {params.finder}/extract_trna_region.py --ffn {input.fasta} --trnas {input.trnas} --out_fa {output.out_fa} > {log} 2>&1
        """

rule blast_mge:
    input:
        fna=os.path.join(GENOMES_DIR, "{sample}.fna"),
        query=os.path.join(RESULTS_DIR, "{sample}", "mge_query.fa"),
        search_params=SEARCH_PARAMS,
    output:
        blast_raw=os.path.join(RESULTS_DIR, "{sample}", "mge_blast_raw.tsv"),
    log:
        os.path.join(RESULTS_DIR, "{sample}", "blast_mge.log")
    params:
        finder=FINDER,
        tmp_dir=os.path.join(RESULTS_DIR, "{sample}"),
    shell:
        """
        python {params.finder}/annotate_mge_region.py --ffn {input.fna} --query {input.query} --out_tsv {output.blast_raw} --tmp_dir {params.tmp_dir} --params {input.search_params} > {log} 2>&1
        """

# Same cascade as the manuscript finder, without the cohort deduplication:
# exact duplication, integrase length, attL non-coding, attL length, and
# optionally no alignment gaps. Thresholds are in search_params.yaml.
rule filter_confident_ie:
    input:
        trna=os.path.join(RESULTS_DIR, "{sample}", "integrase_trna.tsv"),
        integrase_hits=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits_summary.tsv"),
        blast_raw=os.path.join(RESULTS_DIR, "{sample}", "mge_blast_raw.tsv"),
        orfs_gff=os.path.join(RESULTS_DIR, "{sample}", "orfs.gff"),
        fasta=os.path.join(GENOMES_DIR, "{sample}.fna"),
        search_params=SEARCH_PARAMS,
    output:
        audit=os.path.join(RESULTS_DIR, "{sample}", "ie_filter_audit.tsv"),
    log:
        os.path.join(RESULTS_DIR, "{sample}", "filter_confident_ie.log")
    params:
        finder=FINDER
    shell:
        """
        python {params.finder}/filter_confident_ie.py \
            --sample {wildcards.sample} \
            --trna {input.trna} \
            --integrase-hits {input.integrase_hits} \
            --blast-raw {input.blast_raw} \
            --orfs-gff {input.orfs_gff} \
            --fasta {input.fasta} \
            --params {input.search_params} \
            --out-audit {output.audit} \
            > {log} 2>&1
        """

rule export_genome_features:
    input:
        fasta=os.path.join(GENOMES_DIR, "{sample}.fna"),
        audit=os.path.join(RESULTS_DIR, "{sample}", "ie_filter_audit.tsv"),
        integrases=os.path.join(RESULTS_DIR, "{sample}", "integrase_hits_summary.tsv"),
        trna=os.path.join(RESULTS_DIR, "{sample}", "integrase_trna.tsv"),
    output:
        gff=os.path.join(RESULTS_DIR, "{sample}", "{sample}.ie.gff3"),
        gbk=os.path.join(RESULTS_DIR, "{sample}", "{sample}.ie.gbk"),
        sites=os.path.join(RESULTS_DIR, "{sample}", "attachment_sites_genome.tsv"),
    log:
        os.path.join(RESULTS_DIR, "{sample}", "export_genome_features.log")
    params:
        site=os.path.join(workflow.basedir, "scripts"),
        extra="--all-candidates" if ALL_CANDIDATES else "",
    shell:
        """
        python {params.site}/export_genome_features.py \
            --fasta {input.fasta} \
            --audit {input.audit} \
            --integrases {input.integrases} \
            --trna {input.trna} \
            --out-gff3 {output.gff} \
            --out-gbk {output.gbk} \
            --out-tsv {output.sites} \
            --sample {wildcards.sample} \
            {params.extra} \
            > {log} 2>&1
        """
