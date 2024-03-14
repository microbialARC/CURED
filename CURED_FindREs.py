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

PARSER = argparse.ArgumentParser(prog='CURED_FindREs.py', description='This script is a part of the CURED pipeline. This script is used to find restriction enzyme sites in the identified k-mers.')

PARSER.add_argument('--case_control_file',
type=str, help='Csv file of cases and controls.')
PARSER.add_argument('--specificity', '-S',
type=float, help='Specificity for finding RE sites in controls. Default = 100.')
PARSER.add_argument('--added_bases',
type=int, help='Number of added bases on either end of the sequence. Default = 20.')
PARSER.add_argument('--min_coverage',
type=int, help='Minimum coverage threshold for sequence to be considered found in controls. Default = 90.')
PARSER.add_argument('--extension', '-x', type=str,
help='Extension of input assemblies. Default = fna')
PARSER.add_argument('--pcr_product_upstream',
'-UP', type=int, help='Number of bases to include upstream of identified k-mer in outputted PCR product.')
PARSER.add_argument('--pcr_product_downstream',
'-DOWN', type=int, help='Number of bases to include downstream of identified k-mer in outputted PCR product.')
PARSER.add_argument('--compare_coordinates',action='store_true',help='Mode to compare RE by position to determine uniqueness. Default is to determine uniqueness based on presence/absence.')
PARSER.add_argument('--enzymes',type=str,help='Provide a file of restriction enzymes to be used. Default is all enzymes supplied by NE Biolabs.')
PARSER.add_argument('kmers',
type=str, help='List of kmers to be searched.')
PARSER.add_argument('genomes_folder',
type=str, help='Path to genomes.')

ARGS = PARSER.parse_args()
###############################################################################
if ARGS.specificity or ARGS.specificity == 0:
    specificity = ARGS.specificity
else:
    specificity = 100
if ARGS.min_coverage:
    cov_threshold = (ARGS.min_coverage / 100)
else:
    cov_threshold = .90
if ARGS.added_bases or ARGS.added_bases == 0:
    added_bases = ARGS.added_bases
else:
    added_bases = 20
if not ARGS.case_control_file:
    print('Please specifiy case and control file. Exiting script...')
    sys.exit(0)
else:
    accession_numbers = ARGS.case_control_file
if ARGS.pcr_product_upstream:
    include_upstream = ARGS.pcr_product_upstream
else:
    include_upstream = 500
if ARGS.pcr_product_downstream:
    include_downstream = ARGS.pcr_product_downstream
else:
    include_downstream = 500
###############################################################################
TIMESTR = time.strftime("%Y%m%d_%H%M%S")
LOG_FILE = 'CURED_REs_'+ TIMESTR
LOG_LIST = [
    logging.FileHandler("{}.log".format(LOG_FILE)),
    logging.StreamHandler()
]
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s",
    handlers=LOG_LIST)
###############################################################################
# Check dependencies
try:
    blast_check = subprocess.run(['blastn', '-h'], capture_output=True, text=True)
    logging.info('Found blastn.')
except FileNotFoundError:
    logging.error('blast cannot be found. Please install.\n')
    PARSER.exit(status=0)
try:
    samtools_check = subprocess.run(['samtools', '-h'], capture_output=True, text=True)
    logging.info('Found samtools.')
except FileNotFoundError:
    logging.error('samtools cannot be found. Please install.\n')
    PARSER.exit(status=0)
try:
    bwa_check = subprocess.run(['bwa','mem', '-help'], capture_output=True, text=True)
    logging.info('Found bwa.')
except FileNotFoundError:
    logging.error('bwa cannot be found. Please install.\n')
    PARSER.exit(status=0)
try:
    import pkg_resources
    pkg_resources.get_distribution('biopython')
    logging.info("Biopython is installed.")
except pkg_resources.DistributionNotFound:
    logging.error('Biopython cannot be found. Please install.\n')
    PARSER.exit(status=0)

###############################################################################
if ARGS.extension:
    EXT = '.' + ARGS.extension
else:
    EXT = '.fna'

genome_dict = {}
for line in open(accession_numbers, 'r', encoding='utf-8-sig'):
    line = line.rstrip()
    line_list = line.split(',')
    genome = line_list[0]
    classification = line_list[1]
    if classification == 'case':
        genome_dict[genome] = 'CASE'
    else:
        if classification == 'control':
            genome_dict[genome] = 'CONTROL'

value_counter = Counter(genome_dict.values())
num_controls = value_counter['CONTROL']
num_cases = value_counter['CASE']

specificity = specificity / 100
specificity = round(specificity, 2)

# Find maximum nuber of controls that can be found
specificity_step1 = num_controls * specificity # 70000 * .9 = 63,0000

specificity_step2 = int(num_controls - specificity_step1) # 70,000 - 63,000 = 7000

global_specificity = specificity_step2
logging.info(f'Maximum number of controls that a RE can be found in: {str(global_specificity)}')

# This will come from the unitig_master_dict output from the main script
kmer_list = ARGS.kmers
unique_kmers = []
for line in open(kmer_list, 'r'):
    line = line.rstrip()
    kmer = line
    unique_kmers.append(kmer)

# Report to be written to
unique_enzyme_report = open('CURED_UniqueEnzymes.tsv', 'w')
unique_enzyme_report.write('Unique Kmer\tEnzyme Name\tNumber of Controls Enzyme is Found In\tSpecificity\n')

# Report to write PCR product to. (Originally, was supposed to be 400 region)
unique_case_regions = open('CURED_UniqueEnzymes_PCR_Products.tsv', 'w')

# Summary report to write statistics to.
summary_report = open('CURED_FindREs_summary.txt', 'w')
summary_report.write('Kmer\tCase Genome Used\tNumber of Control Genomes Checked\tNumber of Control Genomes Excluded Because Alignment Not Found\tNumber of Control Genomes Excluded Because Minimum Coverage Not Met\n')

# Control file Report
control_report = open('CURED_FindREs_controls.txt', 'w')
control_report.write('Kmer\tEnzyme\tControls Found With Enzyme\n')

# Where the genome files are
extractedFiles_directory = ARGS.genomes_folder

# Enzymes to be searched. Default is what is provided by NEB.
enzyme_query = []
if ARGS.enzymes:
    for line in open(ARGS.enzymes, 'r'):
        line = line.rstrip()
        enzyme_query.append(line)
else:
    enzyme_query = 'N'

# Function to remove temporary files that are created when finding unique restriction sites.
def remove_files_with_pattern(directory, pattern):
   matching_files = glob.glob(os.path.join(directory, pattern))

   for file_path in matching_files:
       try:
           os.remove(file_path)
       except OSError as e:
           logging.warning(f"Error while removing {file_path}: {e}")

# Function to get length of sequence from CIGAR string
def extract_length(cigar_string):
    total_length = 0
    current_length = ""

    for char in cigar_string:
        if char.isdigit():
            current_length += char
        else:
            if current_length:
                total_length += int(current_length)
                current_length = ""

    return total_length

# Function to use Biopython to find REs
def restriction_enzyme_analysis(kmer, enzymes):
    kmer = Seq(kmer)
    try:
        if isinstance(enzymes, set) or isinstance(enzymes, list):
            rb = RestrictionBatch(enzymes)
        else:
            rb = RestrictionBatch(first=[],suppliers='N')
    except Exception as e:
        logging.error('There was an issue with creating the RestrictionBatch. Please ensure that the file is formatted properly.')
        sys.exit(0)

    analysis = Analysis(rb, kmer)
    dict_of_found_enzymes = analysis.with_sites()
    return dict_of_found_enzymes

# Function to add case sequence to PCR dict, to be printed later if sequence has unique RE site.
def add_to_pcr_dict(pcr_dict, kmer, start, end, include_upstream, include_downstream, contig):
    start = start - include_upstream
    end = end + include_downstream
    pcr_dict[kmer] = [start, end, contig]

# Function to run samtools faidx
def samtools_extract(filepath, contig, start, end, rc=None):
    print(rc)
    try:
        sequence_list = []
        if rc == None:
            samtools_extract_cmd = f'samtools faidx {filepath} "{contig}":{start}-{end}'
        else:
            if rc == '-i':
                samtools_extract_cmd = f'samtools faidx {filepath} "{contig}":{start}-{end} -i'
        run_samtools = subprocess.run(samtools_extract_cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        logging.info(samtools_extract_cmd)
        samtools_out = run_samtools.stdout.rstrip().split('\n')
        for line in samtools_out[1:]:
            line = line.rstrip()
            sequence_list.append(line)
        full_sequence = ''.join(sequence_list)
        return full_sequence
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred while running samtools: {e}")
        return None

# Function to parse the output from bwa mem
def parse_bwa_output(bwa_case_output_lines, filepath, unique_kmer):
    try:
        for line in bwa_case_output_lines:
            if not line.startswith('@'):
                line_list = line.split('\t')
                status = line_list[1]
                if status != '4': # In bwa, a flag of 4 means that no alignment was able to be found
                    length_str = line_list[5]
                    extracted_length = extract_length(length_str)
                    end_pos = (extracted_length - 1) + int(line_list[3])
                    end_adjusted = end_pos + added_bases
                    start_adjusted = int(line_list[3]) - added_bases
                    if start_adjusted < 0:
                        start_adjusted = 1
                    aln_length = end_adjusted - start_adjusted
                    contig = line_list[2]

                    sequence = samtools_extract(filepath, contig, start_adjusted, end_adjusted)

                    filepath_basename = filepath.split('/')[-1]
                    if filepath_basename == case_genome_to_check:
                        add_to_pcr_dict(pcr_dict, unique_kmer, start_adjusted, end_adjusted, include_upstream, include_downstream, contig)
                        return sequence, aln_length
                    else:
                        return sequence, aln_length
                else:
                    return None
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        return None

# Function to run bwa mem
def run_bwa(kmer, filepath):
    try:
        input = f">unique_kmer\n{str(kmer)}"
        matching_files = glob.glob(f'{filepath}.a*')
        if not matching_files:
            bwa_index = f'bwa index {filepath}'
            subprocess.run(bwa_index, shell=True)
            logging.info(bwa_index)
        else:
            logging.info('Index files found. Running bwa mem.')
        bwa_mem = f'bwa mem -M -t 16 -v 1 {filepath} -'
        logging.info(bwa_mem)
        run_bwa = subprocess.run(bwa_mem, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, input=input)
        bwa_output_lines = run_bwa.stdout.rstrip().split('\n')
        return bwa_output_lines
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred while running bwa: {e}")
        return None

def parse_blastn_short(blast_output_WS, added_bases, pcr_dict, filepath, unique_kmer, include_upstream, include_downstream):
    strand = blast_output_WS[-1]
    if strand == 'minus':
        start_adjusted = int(blast_output_WS[3]) - added_bases
        end_adjusted = int(blast_output_WS[2]) + added_bases
    else:
        start_adjusted = int(blast_output_WS[2]) - added_bases
        end_adjusted = int(blast_output_WS[3]) + added_bases

    if start_adjusted < 0:
        start_adjusted = 1

    aln_length = end_adjusted - start_adjusted
    contig = blast_output_WS[1]

    pcr_dict[unique_kmer] = [(start_adjusted - include_upstream), (end_adjusted + include_downstream), contig]
    try:
        if strand == 'minus':
            sequence = samtools_extract(filepath, contig, start_adjusted, end_adjusted, rc='-i')
        else:
            sequence = samtools_extract(filepath, contig, start_adjusted, end_adjusted)
        return sequence, aln_length
    except Exception as e:
        logging.error(f"Error occurred while extracting sequence: {e}")
        return None

def run_blastn_short(kmer, filepath):
    try:
        input = f">unique_kmer\n{str(kmer)}"
        blast_cmd = f'blastn -query - -subject {filepath} -word_size 10 -evalue 1e-2 -outfmt "6 qseqid sseqid sstart send sstrand"'
        run_blast = subprocess.run(blast_cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, input=input)
        logging.info(blast_cmd)
        WS_output_lines = run_blast.stdout.rstrip().split('\n')
        blast_output_WS = WS_output_lines[0].split('\t')
        return blast_output_WS
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred while running samtools: {e}")
        return None

# Function to write unique enzyme report
def write_enzyme_report(unique_kmer,enzyme_dict_counter, num_controls, unique_enzyme_report):
    for enzyme, count in enzyme_dict_counter.items():
        final_sensitivity = count/num_controls
        rounded_sensitivity = 1 - round(final_sensitivity, 2)
        sensitivity_as_percent = rounded_sensitivity * 100
        unique_enzyme_report.write(f'{unique_kmer}\t{enzyme}\t{str(count)}\t{str(sensitivity_as_percent)}\n')

# Function to write control report
def write_control_report(unique_kmer, found_controls, control_report):
    for control_enzyme, control_found in found_controls.items():
        control_report.write(f'{unique_kmer}\t{control_enzyme}\t{" ".join(control_found)}\n')

# Function to write PCR report
def write_PCR_report(pcr_dict, unique_kmer, filepath, unique_case_regions):
    start_coordinate = pcr_dict[unique_kmer][0]
    end_coordinate = pcr_dict[unique_kmer][1]
    contig_name = pcr_dict[unique_kmer][-1]
    if start_coordinate < 1:
        start_coordinate = 1
    extract_pcr_cmd = f'samtools faidx {filepath} "{contig_name}":{start_coordinate}-{end_coordinate} > {unique_kmer}_pcr.fa'
    run_extract_pcr = subprocess.run(extract_pcr_cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    with open(f'{unique_kmer}_pcr.fa', 'r') as pcr_file:
        for line in pcr_file:
            line = line.rstrip()
            if line.startswith('>'):
                header = f'>{unique_kmer}'
                unique_case_regions.write(header + '\n')
            else:
                unique_case_regions.write(line + '\n')
    os.remove(f'{unique_kmer}_pcr.fa')

# Get paths of controls to be used during the for loop used for finding unique restriction sites.
CONTROLS = [genome for genome, status in genome_dict.items() if status == 'CONTROL']

# Get number of controls and the filepaths for each control
controls = len(CONTROLS)
control_filepaths = [os.path.join(extractedFiles_directory, file) for file in CONTROLS if os.path.exists(os.path.join(extractedFiles_directory, file))]

# Get case genome to check
case_genome_to_check = None
for accession, label in genome_dict.items():
    if label == 'CASE':
        case_genome_to_check = accession
        break
logging.info(f'Case genome being used for finding unique restriction enzyme sites: {case_genome_to_check}')

pcr_dict = {} # Initiate PCR dictionary where sequences will be stored

# Enter loop to parse each k-mer one by one.
for unique_kmer in unique_kmers:

    logging.info(f'Kmer being searched: {str(unique_kmer)}')

    dict_of_found_enzymes = restriction_enzyme_analysis(unique_kmer, enzyme_query)

    logging.info(f'Checking for RES in the case genome. Genome being checked: {case_genome_to_check}')
    logging.info(f'Enzymes found in {case_genome_to_check}: {str(dict_of_found_enzymes)}')

    if dict_of_found_enzymes: # If REs were found in the unique kmer, execute this block of code.
        logging.info(f'Matches to an RES in {unique_kmer} were found in {case_genome_to_check}')

        case_enzyme_set = set(dict_of_found_enzymes.keys())
        enzyme_dict_counter = dict.fromkeys(dict_of_found_enzymes.keys(), 0)

        logging.info(f'Initial Enzyme Dictionary: {str(enzyme_dict_counter)}')

        filepath = os.path.join(extractedFiles_directory, case_genome_to_check)
        bwa_output_case = run_bwa(unique_kmer, filepath)
        bwa_output_case_list = parse_bwa_output(bwa_output_case, filepath, unique_kmer)

        if bwa_output_case_list == None:
            blastn_output = run_blastn_short(unique_kmer, filepath)
            case_sequence, aln_length = parse_blastn_short(blastn_output, added_bases, pcr_dict, filepath, unique_kmer, include_upstream, include_downstream)

        else:
            case_sequence = bwa_output_case_list[0]
            aln_length = bwa_output_case_list[1]

        # Calculate minimum coverage to be considered
        min_coverage = int(cov_threshold * int(aln_length))

        no_aln_counter = 0
        min_cov_counter = 0
        num_controls_checked = 0

        found_controls = {}
        for control in control_filepaths:
            control_basename = control.split('/')[-1].split(EXT)[0]

            bwa_output_control = run_bwa(case_sequence, control)
            num_controls_checked += 1

            bwa_output_control_list = parse_bwa_output(bwa_output_control, control, unique_kmer)

            if bwa_output_control_list == None:
                no_aln_counter += 1
                continue
            else:
                control_sequence = bwa_output_control_list[0]
                control_aln_length = bwa_output_control_list[1]

            if int(control_aln_length) <= min_coverage:
                logging.info('Length of alignment in control: ' + str(control_aln_length))
                logging.warning('Alignment in control is low coverage. Searching next...')
                min_cov_counter += 1
                continue

            control_enzymes_found = restriction_enzyme_analysis(control_sequence, case_enzyme_set)
            logging.info(f'Enzymes found in control: {str(control_enzymes_found)}')

            to_be_removed = []

            if ARGS.compare_coordinates:
                for control_E,control_V in control_enzymes_found.items():
                    if control_E in enzyme_dict_counter:
                        for coord in control_V:
                            start_flanking_region = len(control_sequence) - added_bases
                            if (1 <= coord < added_bases) or (start_flanking_region < coord <= len(fcontrol_sequence)):
                                continue
                            else:
                                enzyme_dict_counter[control_E] += 1
                                if control_E in found_controls:
                                    found_controls[control_E].append(control_basename)
                                else:
                                    found_controls[control_E] = [control_basename]
                                break
                    else:
                        continue

                for E,counter in enzyme_dict_counter.items():
                    if counter > specificity_step2:
                            to_be_removed.append(E)

            if not ARGS.compare_coordinates:
                for enzyme in control_enzymes_found.keys():
                    if enzyme in enzyme_dict_counter:
                        enzyme_dict_counter[enzyme] += 1
                        if enzyme in found_controls:
                            found_controls[enzyme].append(control_basename)
                        else:
                            found_controls[enzyme] = [control_basename]

                for E,counter in enzyme_dict_counter.items():
                    if counter > specificity_step2:
                        to_be_removed.append(E)

            logging.info(f'Enzyme Counter Dictionary: {str(enzyme_dict_counter)}')

            for element in to_be_removed:
                enzyme_dict_counter.pop(element)

            case_enzyme_set.difference_update(set(to_be_removed))

            logging.info(f'Updated enzyme counter dictionary: {str(enzyme_dict_counter)}')

            if len(case_enzyme_set) == 0:
                logging.info('Reached maximum number of controls that enzyme can be found in. Searching next...')
                break
            else:
                continue

        # Write to Enzyme Report
        write_enzyme_report(unique_kmer,enzyme_dict_counter, num_controls, unique_enzyme_report)

        # Write to Control Report
        write_control_report(unique_kmer, found_controls, control_report)

        # Write to PCR Product Report
        if enzyme_dict_counter:
            write_PCR_report(pcr_dict, unique_kmer, filepath, unique_case_regions)

        summary_report.write(f'{unique_kmer}\t{case_genome_to_check}\t{num_controls_checked}\t{no_aln_counter}\t{min_cov_counter}\n')

        logging.info('Deleting tmp files. Moving to next kmer.')
        logging.info('-----------------------------------------------')

    else: # If not enzyme_occurrences
        logging.info('No Restriction Sites found in this unique Kmer.')
        logging.info('-----------------------------------------------')
        continue
