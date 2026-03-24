#!/usr/bin/env python3
import re
import os
import sys
import subprocess
import glob
import time
from Bio.Restriction import *
from Bio.Seq import Seq
from collections import Counter
import argparse
import logging
import hashlib
from Bio import SeqIO

# =============================================================================
# CONFIGURATION & ARGUMENTS
# =============================================================================
PARSER = argparse.ArgumentParser(prog='CURED_FindREs_v25.py', description='Find unique RE sites using Direct Extraction & Unified Minimap2.')

PARSER.add_argument('--case_control_file', type=str, required=True, help='Csv file of cases and controls.')
PARSER.add_argument('--specificity', '-S', type=float, help='Specificity for finding RE sites in controls. Default = 100.')
PARSER.add_argument('--added_bases', type=int, help='Number of added bases on either end of the sequence. Default = 20.')
PARSER.add_argument('--min_coverage', type=int, help='Minimum coverage threshold (0-100). Default = 90.')
PARSER.add_argument('--extension', '-x', type=str, help='Extension of input assemblies. Default = fna')
PARSER.add_argument('--pcr_product_upstream', '-UP', type=int, help='Upstream bases for PCR product.')
PARSER.add_argument('--pcr_product_downstream', '-DOWN', type=int, help='Downstream bases for PCR product.')
PARSER.add_argument('--compare_coordinates', action='store_true', help='Compare RE by position.')
PARSER.add_argument('--tolerance', type=int, default=15, help='Resolution tolerance (bp) for compare_coordinates. Default = 15.')
PARSER.add_argument('--output_dir', '-o', type=str, default='CURED_Results', help='Directory to save all output files.')
PARSER.add_argument('--enzymes', type=str, help='File of restriction enzymes. Default: All NEB.')
PARSER.add_argument('--similarity_matrix', type=str, help='Path to similarity matrix csv.')
PARSER.add_argument('--batch_size', '-b', type=int, default=5000, help='Batch size. Default = 5000.')
PARSER.add_argument('--threads', '-t', type=int, default=38, help='Threads. Default = 38.')
PARSER.add_argument('kmers', type=str, help='List of kmers to be searched.')
PARSER.add_argument('genomes_folder', type=str, help='Path to genomes.')

# MINIMAP2 PARAMETERS
PARSER.add_argument('--k_size', '-k', type=int, default=11, help='Minimap2 seed size. Default = 11.')
PARSER.add_argument('--mm2_m', '-m', type=int, default=10, help='Minimap2 -m. Default=10.')
PARSER.add_argument('--mm2_s', '-s', type=int, default=10, help='Minimap2 -s. Default=10.')
PARSER.add_argument('--mm2_A', '-A', type=int, default=2, help='Minimap2 -A. Default=2.')
PARSER.add_argument('--mm2_B', '-B', type=int, default=4, help='Minimap2 -B. Default=4.')
PARSER.add_argument('--mm2_w', '-w', type=int, default=10, help='Minimap2 -w. Default=10.')

ARGS = PARSER.parse_args()

# CREATE OUTPUT DIRECTORY
if not os.path.exists(ARGS.output_dir):
    os.makedirs(ARGS.output_dir)
    print(f"Created output directory: {ARGS.output_dir}")
else:
    print(f"Using existing output directory: {ARGS.output_dir}")

# CONSTRUCT COMMANDS
SAFE_FLAGS = "-N 200 -p 0.01 -f 5000"
MM2_PARAMS_SAFE = (f"-c -x sr -k {ARGS.k_size} -w {ARGS.mm2_w} -m {ARGS.mm2_m} "
                   f"-s {ARGS.mm2_s} -A {ARGS.mm2_A} -B {ARGS.mm2_B} {SAFE_FLAGS} --cs")
#MM2_PARAMS_CASE = (f"-c -x sr -k {ARGS.k_size} -w {ARGS.mm2_w} -m {ARGS.mm2_m} "
#                   f"-s {ARGS.mm2_s} -A {ARGS.mm2_A} -B {ARGS.mm2_B} -N 1 --cs")
MM2_PARAMS_CASE = "-c -x sr -k 7 -w 1 -m 10 -s 10 -A 2 -B 4 -N 1 --cs"
# Setup Globals
specificity = ARGS.specificity if ARGS.specificity is not None else 100
cov_threshold = (ARGS.min_coverage / 100) if ARGS.min_coverage else 0.90
added_bases = ARGS.added_bases if ARGS.added_bases is not None else 20
batch_size = ARGS.batch_size if ARGS.batch_size else 5000
num_threads = ARGS.threads if ARGS.threads else 38
tolerance_bp = ARGS.tolerance

accession_numbers = ARGS.case_control_file
    
include_upstream = ARGS.pcr_product_upstream if ARGS.pcr_product_upstream else 500
include_downstream = ARGS.pcr_product_downstream if ARGS.pcr_product_downstream else 500

TIMESTR = time.strftime("%Y%m%d_%H%M%S")
LOG_FILE = os.path.join(ARGS.output_dir, 'CURED_REs_'+ TIMESTR + '.log')
logging.basicConfig(level=logging.INFO, format="%(asctime)s:%(levelname)s:%(message)s", handlers=[
    logging.FileHandler(LOG_FILE),
    logging.StreamHandler()
])

logging.info("="*60)
logging.info(f"COMMAND EXECUTED:  python {' '.join(sys.argv)}")
logging.info(f"WORKING DIRECTORY: {os.getcwd()}")
logging.info(f"CONDA ENVIRONMENT: {os.environ.get('CONDA_DEFAULT_ENV', 'Not Detected')}")
logging.info(f"PARSED ARGUMENTS:  {ARGS}")
logging.info("="*60)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def get_batches(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def check_dep(cmd, name):
    try:
        subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        logging.error(f'{name} cannot be found. Please install.')
        sys.exit(0)

check_dep(['minimap2', '--version'], 'minimap2')
check_dep(['samtools', '--version'], 'samtools')

def samtools_extract(filepath, contig, start, end, rc=None):
    try:
        sequence_list = []
        cmd = f'samtools faidx {filepath} "{contig}":{start}-{end}'
        if rc == '-i': cmd += ' -i'
        run = subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        lines = run.stdout.rstrip().split('\n')
        for line in lines[1:]:
            sequence_list.append(line.rstrip())
        return ''.join(sequence_list)
    except subprocess.CalledProcessError as e:
        return None

def restriction_enzyme_analysis(kmer, enzymes):
    kmer = Seq(kmer)
    try:
        if isinstance(enzymes, list): rb = RestrictionBatch(enzymes)
        else: rb = RestrictionBatch(first=[],suppliers='N')
    except: return {}
    return Analysis(rb, kmer).with_sites()

def write_enzyme_report(unique_kmer, enzyme_dict_counter, num_controls, unique_enzyme_report):
    for enzyme, count in enzyme_dict_counter.items():
        sens = (1 - (count/num_controls)) * 100
        unique_enzyme_report.write(f'{unique_kmer}\t{enzyme}\t{str(count)}\t{str(sens)}\n')

def write_control_report(unique_kmer, found_controls, control_report):
    for control_enzyme, control_found in found_controls.items():
        control_report.write(f'{unique_kmer}\t{control_enzyme}\t{" ".join(control_found)}\n')

# CRASH-PROOF PCR REPORT
def write_PCR_report(pcr_dict, unique_kmer, filepath, unique_case_regions):
    if unique_kmer not in pcr_dict: return
    start, end, contig = pcr_dict[unique_kmer]
    start = max(1, start)
    khash = hashlib.md5(unique_kmer.encode()).hexdigest()
    
    # Save temp PCR file inside output dir to stay clean
    temp_pcr_path = os.path.join(ARGS.output_dir, f'{khash}_pcr.fa')
    
    # 1. Ensure Index Exists (Silent Fix)
    if not os.path.exists(filepath + ".fai"):
        subprocess.run(f'samtools faidx {filepath}', shell=True)

    # 2. Run Extraction with CHECK
    cmd = f'samtools faidx {filepath} "{contig}":{start}-{end} > {temp_pcr_path}'
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError:
        logging.warning(f"Failed to generate PCR product for {unique_kmer}. Skipping.")
        return

    # 3. Read File SAFELY (No Crash)
    if os.path.exists(temp_pcr_path):
        if os.path.getsize(temp_pcr_path) > 0:
            with open(temp_pcr_path, 'r') as f:
                for line in f:
                    if line.startswith('>'): unique_case_regions.write(f'>{unique_kmer}\n')
                    else: unique_case_regions.write(line)
        os.remove(temp_pcr_path)
    else:
        logging.warning(f"PCR File missing for {unique_kmer} even after Samtools run.")

# =============================================================================
# MAIN SCRIPT EXECUTION
# =============================================================================
if ARGS.extension: EXT = '.' + ARGS.extension
else: EXT = '.fna'

genome_dict = {}
for line in open(accession_numbers, 'r', encoding='utf-8-sig'):
    line = line.rstrip()
    line_list = line.split(',')
    genome = line_list[0]
    if not genome.endswith(EXT): genome = genome + EXT
    if line_list[1] == 'case': genome_dict[genome] = 'CASE'
    else: genome_dict[genome] = 'CONTROL'

value_counter = Counter(genome_dict.values())
num_controls = value_counter['CONTROL']
specificity_adj = round(specificity / 100, 2)
global_specificity = int(num_controls * (1 - specificity_adj))
logging.info(f'Max Allowed Control Hits: {global_specificity}')

unique_kmers = []
for line in open(ARGS.kmers, 'r'):
    unique_kmers.append(line.strip())

# File Handles (ALL IN OUTPUT DIR)
path_unique_enzyme = os.path.join(ARGS.output_dir, 'CURED_UniqueEnzymes_minimap2.tsv')
path_unique_case_regions = os.path.join(ARGS.output_dir, 'CURED_UniqueEnzymes_PCR_Products.tsv')
path_summary_report = os.path.join(ARGS.output_dir, 'CURED_FindREs_summary.txt')
path_control_report = os.path.join(ARGS.output_dir, 'CURED_FindREs_controls.txt')

unique_enzyme_report = open(path_unique_enzyme, 'w')
unique_enzyme_report.write('Unique Kmer\tEnzyme Name\tNumber of Controls Enzyme is Found In\tSpecificity\n')

unique_case_regions = open(path_unique_case_regions, 'w')

summary_report = open(path_summary_report, 'w')
summary_report.write('Kmer\tCase Genome Used\tNumber of Control Genomes Checked\tStats\n')

control_report = open(path_control_report, 'w')
control_report.write('Kmer\tEnzyme\tControls Found With Enzyme\n')

extractedFiles_directory = ARGS.genomes_folder
enzyme_query = []
if ARGS.enzymes:
    for line in open(ARGS.enzymes, 'r'):
        enzyme_query.append(line.strip())
else:
    enzyme_query = 'N'

# --- SIMILARITY MATRIX LOGIC ---
if ARGS.similarity_matrix:
    sim_mat = {}
    for line in open(ARGS.similarity_matrix, 'r', encoding='utf-8-sig'):
        line = line.rstrip()
        line_list = line.split(',')
        genome = line_list[0]
        sim_count = line_list[1]
        if not genome.endswith(EXT):
            genome = genome + EXT
        sim_mat[genome] = sim_count

    num_genomes = len(sim_mat)
    midpoint = num_genomes // 2
    similarity_tracked_genomes = list(sim_mat.keys())
    if num_genomes > 200:
        CONTROLS = similarity_tracked_genomes[:50] + similarity_tracked_genomes[midpoint-25:midpoint+25] + similarity_tracked_genomes[-50:] + similarity_tracked_genomes[50:midpoint-25] + similarity_tracked_genomes[midpoint+25:-50]
    else:
        CONTROLS = [genome for genome, status in genome_dict.items() if status == 'CONTROL']
else:
    CONTROLS = [genome for genome, status in genome_dict.items() if status == 'CONTROL']
# --------------------------------

control_filepaths = [os.path.join(extractedFiles_directory, f) for f in CONTROLS if os.path.exists(os.path.join(extractedFiles_directory, f))]
case_genome_to_check = next((g for g, s in genome_dict.items() if s == 'CASE'), None)
logging.info(f'Case genome: {case_genome_to_check}')

# ------------------------------------------------------------------------------
# PHASE 1: PRE-PROCESSING (Minimap2 Mapping to Case)
# ------------------------------------------------------------------------------
logging.info("PHASE 1: Mapping all k-mers to Case Genome (Using Minimap2)...")
pcr_dict = {}
query_data = {}
kmers_orig = []
queries_fasta_path = os.path.join(ARGS.output_dir, "all_queries.fasta") 

raw_kmers_fasta = os.path.join(ARGS.output_dir, "raw_kmers_temp.fasta") 
with open(raw_kmers_fasta, 'w') as f:
    for k in unique_kmers:
        f.write(f">{k}\n{k}\n")
        kmers_orig.append(k)

case_path_full = os.path.join(extractedFiles_directory, case_genome_to_check)
mm2_case_cmd = f"minimap2 -t {num_threads} {MM2_PARAMS_CASE} {case_path_full} {raw_kmers_fasta}"

logging.info(f"  - Cmd: {mm2_case_cmd}")
proc_case = subprocess.Popen(mm2_case_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

mapped_count = 0
with open(queries_fasta_path, 'w') as q_out:
    for line in proc_case.stdout:
        cols = line.strip().split('\t')
        if len(cols) < 10: continue
        
        kmer_name = cols[0]
        contig = cols[5]
        t_start = int(cols[7])
        t_end = int(cols[8])
        
        start_adj = max(1, t_start - added_bases)
        end_adj = t_end + added_bases
        
        extended_seq = samtools_extract(case_path_full, contig, start_adj, end_adj)
        if not extended_seq: continue

        if kmer_name not in kmers_orig:
            kmer_name = str(Seq(kmer_name).reverse_complement())
        
        initial_enzymes_dict = restriction_enzyme_analysis(kmer_name, enzyme_query)
        if not initial_enzymes_dict: continue
        
        # --- OPTIMIZATION: PRE-COMPILE RestrictionBatch ---
        # Ensure keys are STRINGS
        relevant_enz_names = [str(enz) for enz in initial_enzymes_dict.keys()]
        
        try:
            if relevant_enz_names:
                optimized_rb = RestrictionBatch(relevant_enz_names)
            else:
                optimized_rb = None
        except: optimized_rb = None
        # --------------------------------------------------

        # Compare Coordinates Setup
        # FIX IN V25: Use simpler iteration to guarantee population
        enzyme_locs = {}
        if ARGS.compare_coordinates:
            try:
                ext_seq_obj = Seq(extended_seq)
                if optimized_rb:
                    locs_all = optimized_rb.search(ext_seq_obj)
                    # Iterate the RESULTS, not the names, to avoid getattr issues
                    for enz_obj, sites in locs_all.items():
                        enzyme_locs[str(enz_obj)] = sites
            except: pass

        query_data[kmer_name] = {
            'case_seq': extended_seq,
            'case_enzymes': set(relevant_enz_names), 
            'case_enzyme_locs': enzyme_locs,
            'control_hits': {enz: 0 for enz in relevant_enz_names}, 
            'control_names': {enz: set() for enz in relevant_enz_names}, 
            'total_genome_hits': 0,
            'optimized_rb': optimized_rb
        }
        
        pcr_dict[kmer_name] = [start_adj - include_upstream, end_adj + include_downstream, contig]
        
        q_out.write(f">{kmer_name}\n{extended_seq}\n")
        mapped_count += 1

proc_case.wait()
if os.path.exists(raw_kmers_fasta): os.remove(raw_kmers_fasta)
logging.info(f"  - Mapped {mapped_count} K-mers to Case Genome.")

# --- MISSING K-MER CHECK ---
all_input_kmers = set(unique_kmers)
found_kmers = set(query_data.keys())
missing_kmers = all_input_kmers - found_kmers

if missing_kmers:
    print("\n" + "!"*60)
    print("WARNING: Some K-mers were NOT found in the Case Genome!")
    print(f"Input: {len(unique_kmers)} -> Found: {len(found_kmers)}")
    print("The following K-mers will be skipped:")
    for mk in list(missing_kmers)[:10]:
        print(f"  - {mk}")
    if len(missing_kmers) > 10:
        print(f"  ... and {len(missing_kmers)-10} more.")
    print("!"*60 + "\n")
    logging.warning(f"Missing {len(missing_kmers)} K-mers from Case Genome.")
else:
    logging.info("SUCCESS: All input K-mers were found in the Case Genome.")
# --------------------------------

# ------------------------------------------------------------------------------
# PHASE 2: BATCH SCANNING (Direct Extraction - OPTIMIZED)
# ------------------------------------------------------------------------------
BATCH_SIZE = batch_size
logging.info(f"PHASE 2: Scanning {len(control_filepaths)} controls in batches of {BATCH_SIZE}...")
temp_query = os.path.join(ARGS.output_dir, "temp_batch_controls.fasta") 
total_controls_checked = 0

for batch in get_batches(control_filepaths, BATCH_SIZE):
    seq_db = {}
    with open(temp_query, 'w') as t_out:
        for c_file in batch:
            c_name = os.path.basename(c_file).replace(EXT, '')
            try:
                for record in SeqIO.parse(c_file, "fasta"):
                    unique_id = f"{c_name}|{record.id}"
                    t_out.write(f">{unique_id}\n{str(record.seq)}\n")
                    seq_db[unique_id] = str(record.seq)
            except Exception: pass

    mm2_cmd = f"minimap2 -t {num_threads} {MM2_PARAMS_SAFE} {queries_fasta_path} {temp_query}"
    
    process = subprocess.Popen(mm2_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    for line in process.stdout:
        cols = line.strip().split('\t')
        if len(cols) < 10: continue
        
        query_header = cols[0]
        kmer_name = cols[5]
        
        if kmer_name not in query_data: continue
        
        q_start = int(cols[2])
        q_end = int(cols[3])
        strand = cols[4]
        
        if query_header in seq_db:
            full_contig = seq_db[query_header]
            extracted_seq = full_contig[q_start:q_end]
            
            if strand == '-':
                control_seq = str(Seq(extracted_seq).reverse_complement())
            else:
                control_seq = extracted_seq
            
            expected_len = len(query_data[kmer_name]['case_seq'])
            if len(control_seq) < (expected_len * cov_threshold): continue
            
            query_data[kmer_name]['total_genome_hits'] += 1
            
            # --- OPTIMIZED ENZYME ANALYSIS ---
            rb = query_data[kmer_name]['optimized_rb']
            if not rb: continue
            
            try:
                control_seq_obj = Seq(control_seq)
                analysis = rb.search(control_seq_obj)
                
                for enz, found_locs in analysis.items():
                    enz_name = str(enz)
                    if not found_locs: continue
                    
                    if ARGS.compare_coordinates:
                        expected_locs = query_data[kmer_name]['case_enzyme_locs'].get(enz_name, [])
                        found_valid = False
                        for loc in found_locs:
                            if any(abs(loc - exp) <= tolerance_bp for exp in expected_locs):
                                found_valid = True
                                break
                        if not found_valid: continue 

                    query_data[kmer_name]['control_hits'][enz_name] += 1
                    control_genome_name = query_header.split('|')[0]
                    
                    query_data[kmer_name]['control_names'][enz_name].add(control_genome_name)
                    
                    if query_data[kmer_name]['control_hits'][enz_name] > global_specificity:
                         if enz_name in query_data[kmer_name]['case_enzymes']:
                             query_data[kmer_name]['case_enzymes'].remove(enz_name)
            except: pass
            # ---------------------------------

    process.wait()
    if process.returncode != 0:
        error_msg = process.stderr.read() if process.stderr else "No stderr captured."
        print(f"\nCRITICAL ERROR: Minimap2 crashed on batch {total_controls_checked}-{total_controls_checked+len(batch)}.")
        sys.exit(1)

    total_controls_checked += len(batch)
    logging.info(f"    Processed {total_controls_checked} / {len(control_filepaths)} genomes...")
    if os.path.exists(temp_query): os.remove(temp_query)

if os.path.exists(queries_fasta_path): os.remove(queries_fasta_path)

# ------------------------------------------------------------------------------
# PHASE 3: REPORTING
# ------------------------------------------------------------------------------
logging.info("PHASE 3: Generating reports...")

summary_path = os.path.join(ARGS.output_dir, "FINAL_SUMMARY.csv")
stats_total_kmers_input = len(unique_kmers)
stats_kmers_found_in_case = len(query_data)
stats_kmers_with_valid_enzymes = 0
stats_total_kmer_enzyme_pairs = 0

with open(summary_path, "w") as f:
    f.write("Kmer_Name,Total_Enzymes_Case,Unique_Enzymes_Count,Unique_Enzymes_List,Control_Genomes_Hit\n")
    for kmer, data in query_data.items():
        # Force str(enz)
        final_valid_enzymes = [str(enz) for enz, count in data['control_hits'].items() if count <= global_specificity]
        
        if final_valid_enzymes:
            stats_kmers_with_valid_enzymes += 1
            stats_total_kmer_enzyme_pairs += len(final_valid_enzymes)

        total_enz = len(data['control_hits'])
        unique_cnt = len(final_valid_enzymes)
        unique_list_str = ";".join(final_valid_enzymes)
        f.write(f"{kmer},{total_enz},{unique_cnt},{unique_list_str},{data['total_genome_hits']}\n")

for kmer, data in query_data.items():
    # Force str(enz)
    final_enzymes = {str(enz): c for enz, c in data['control_hits'].items() if c <= global_specificity}
    
    write_enzyme_report(kmer, final_enzymes, num_controls, unique_enzyme_report)
    
    found_controls_map = {enz: data['control_names'][enz] for enz in final_enzymes}
    write_control_report(kmer, found_controls_map, control_report)
    
    if final_enzymes:
        filepath = os.path.join(extractedFiles_directory, case_genome_to_check)
        write_PCR_report(pcr_dict, kmer, filepath, unique_case_regions)

    summary_report.write(f'{kmer}\t{case_genome_to_check}\t{total_controls_checked}\tN/A\n')

unique_enzyme_report.close()
unique_case_regions.close()
summary_report.close()
control_report.close()

print("\n" + "="*60)
print("             CURED: FINAL SUMMARY REPORT")
print("="*60)
print(f"Total K-mers Input:              {stats_total_kmers_input}")
print(f"Found in Case Genome:            {stats_kmers_found_in_case}")
print(f"K-mers with Valid Enzymes:       {stats_kmers_with_valid_enzymes}")
print(f"Total Unique Kmer-Enzyme Pairs:  {stats_total_kmer_enzyme_pairs}")
if ARGS.compare_coordinates:
    print(f"Resolution Tolerance Used:       +/- {tolerance_bp} bp")
print(f"Output Directory:                {ARGS.output_dir}")
print(f"CSV Summary saved to:            {summary_path}")
print("="*60 + "\n")

if ARGS.compare_coordinates:
    logging.info(f"Summary: Input={stats_total_kmers_input}, FoundInCase={stats_kmers_found_in_case}, ValidKmers={stats_kmers_with_valid_enzymes}, Pairs={stats_total_kmer_enzyme_pairs}, Tol={tolerance_bp}")
else:
    logging.info(f"Summary: Input={stats_total_kmers_input}, FoundInCase={stats_kmers_found_in_case}, ValidKmers={stats_kmers_with_valid_enzymes}, ValidPairs={stats_total_kmer_enzyme_pairs}")
logging.info("Done.")