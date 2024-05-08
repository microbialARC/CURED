#!/usr/bin/env

import sys
from collections import Counter
import time
import os
import subprocess
import logging
import argparse
import tempfile
import shutil
import zipfile
import json

PARSER = argparse.ArgumentParser(prog='CURED_Main.py',description='This script is part of the CURED pipeline. This script is used for finding unique clonal biomarkers in your cases.')

PARSER.add_argument('--number_of_cases', '-cases',
type=int, help='Add in the number of cases to be used in each iteration.')
PARSER.add_argument('--number_of_controls', '-controls',
type=int, help='Add in the number of controls to be used in each iteration.')
PARSER.add_argument('--threads', '-T',
type=int, help='Number of threads to be used for unitig-caller. Default = 1.')
PARSER.add_argument('--case_control_file',
type=str, help='Csv file of genomes to be used with case or control designation. ')
PARSER.add_argument('--sensitivity',
type=int, help='Specifiy sensitivity. Default = 100.')
PARSER.add_argument('--specificity',
type=int, help='Specify specificity. Default = 100.')
PARSER.add_argument('--kmer_length', '-K',
type=int, help='Specify minimum length of k-mer to search for. Default is 20.')
PARSER.add_argument("--species",
type=str,nargs= '+',help='Genus or species of interest')
PARSER.add_argument("--sequence_type", "-st",
type=int,help='sequence type of interest')
PARSER.add_argument('--case_accession_list', type=str, help='List of case accessions.')
PARSER.add_argument("--extension", "-x",
type=str, help='extension of assembly inputs. Ignore if using --species/--sequence_type options. Default = fna')
PARSER.add_argument('--database', '-db',type=str,help='Choose to download genomes from RefSeq, GenBank, or both. Default = both.')
PARSER.add_argument('--summary', action='store_true', help='Check to see how many genomes will be downloaded if you use the --species option. Use this option with --database and --species options.')
PARSER.add_argument('--quiet', '-q', action='store_true',help='No screen output. Default = OFF')
PARSER.add_argument('--genomes_folder', type=str, help='Path to genomes.')
PARSER.add_argument('--case_genomes', action='store_true',help='Option if you have local sequencing data to serve as the cases. Use --species to download control genomes.')
PARSER.add_argument('--use_datasets', action='store_true',help='Option to provide CURED with a case and control file of ncbi accessions and downloaded them. Use with --case_control_file')
PARSER.add_argument('--use_simple', action='store_true',help='Option to run unitig-caller simple mode. This is useful for when you already have a list of k-mers that you want to query against a set of genomes. ')
PARSER.add_argument('--kmer_list', type=str, help='List of k-mers to be used as query in running unitig-caller simple mode')
PARSER.add_argument('-help', action='help', help='Show this help message and exit.')
ARGS = PARSER.parse_args()

if bool(vars(ARGS)["sequence_type"]) and not bool(vars(ARGS)["species"]):
    PARSER.exit(status=0, message="Error: You have to use --sequence_type with --species\n")
if bool(vars(ARGS)['database']) and not bool(vars(ARGS)['species']):
    PARSER.exit(status=0, message='Error: When using --database option, must use --species\n')
if bool(vars(ARGS)['case_accession_list']) and not bool(vars(ARGS)['species']):
    PARSER.exit(status=0, message='Error: When using --case_accession_list option, must use --species\n')
if bool(vars(ARGS)['species']) and bool(vars(ARGS)['genomes_folder']) and not bool(vars(ARGS)['case_genomes']):
    PARSER.exit(status=0, message='Error: When using --case-genomes option, must use --species and --genomes_folder\n')
if bool(vars(ARGS)['use_datasets']) and not bool(vars(ARGS)['case_control_file']):
    PARSER.exit(status=0, message='Error: When using --use_datasets option, must provide a case/control file using --case_control_file \n')
if bool(vars(ARGS)['use_simple']) and not bool(vars(ARGS)['kmer_list']):
    PARSER.exit(status=0, message='Error: When using --use_simple option, must provide a file of k-mers using --kmer_list\n')
###############################################################################
if ARGS.case_control_file:
    case_control = ARGS.case_control_file
if ARGS.threads:
    thr = ARGS.threads
else:
    thr = 1
if ARGS.kmer_length:
    kmer_length = ARGS.kmer_length
else:
    kmer_length = 20
if ARGS.extension:
    EXT = '.' + ARGS.extension
else:
    EXT = '.fna'
if ARGS.database:
    db = str(ARGS.database)
else:
    db = 'all'
if ARGS.genomes_folder:
    extractedFiles_directory = ARGS.genomes_folder
else:
    extractedFiles_directory = 'genomes/'
if ARGS.species:
    genus_species = ' '.join(ARGS.species)
    genus_species_command = '"' + genus_species + '"'
if ARGS.kmer_list:
    kmer_list = ARGS.kmer_list
if ARGS.sequence_type:
    sequence_type = str(ARGS.sequence_type)

TIMESTR = time.strftime("%Y%m%d_%H%M%S")
LOG_FILE = 'CURED_'+ TIMESTR
if ARGS.quiet:
    LOG_LIST = [
        logging.FileHandler("{}.log".format(LOG_FILE))
    ]
else:
    LOG_LIST = [
        logging.FileHandler("{}.log".format(LOG_FILE)),
        logging.StreamHandler()
    ]
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s",
    handlers=LOG_LIST)

current_directory = os.getcwd()
temp_dir= tempfile.mkdtemp(dir=current_directory)
###############################################################################
# Check dependencies
try:
    unitigcaller_check = subprocess.run(['unitig-caller', '--help'], capture_output=True, text=True)
    logging.info("Found unitig-caller.")
except FileNotFoundError:
    logging.error('ERROR: unitig-caller cannot be found. Please install.')
    PARSER.exit(status=0)
try:
  mlst_version = subprocess.run(['mlst', '-version'], capture_output=True, text=True)
  logging.info(f"Found mlst: {mlst_version.stdout.strip()}")
except FileNotFoundError:
  logging.error('ERROR: mlst cannot be found. Please install.')
  PARSER.exit(status=0)
try:
  datasets_version = subprocess.run(['datasets', '--version'], capture_output=True, text=True)
  logging.info(f"Found ncbi-datasets: {datasets_version.stdout.strip()}")
except FileNotFoundError:
  logging.error('ERROR: ncbi-datasets cannot be found. Please install.')
  PARSER.exit(status=0)
###############################################################################
def extract_files_from_zip(zip_path, file_extension, destination_dir):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        for member in zip_ref.infolist():
            if not member.is_dir() and member.filename.endswith(file_extension):
                filename = os.path.basename(member.filename)
                filename_list = filename.split('_')
                new_filename = '_'.join(filename_list[0:2]) + '.fna' # Check these lines of code
                extracted_path = zip_ref.extract(member, path=destination_dir)
                new_path = os.path.join(destination_dir, new_filename)
                shutil.move(extracted_path, new_path)
                logging.info(f"Moved: {new_path}")

def download_genomes_summary(genus_species_command,db):
    datasets_summary_cmd = subprocess.run([f"datasets summary genome taxon {genus_species_command} --assembly-source {db}"], shell=True, capture_output=True, text=True)
    json_output = datasets_summary_cmd.stdout

    if datasets_summary_cmd.returncode != 0:
            logging.error("Error occurred while executing the command.")
            sys.exit(1)

    data = json.loads(json_output)

    keys = list(data.keys())
    second_key = keys[1] if len(keys) > 1 else None

    second_value = data.get(second_key) if second_key else None
    logging.info(f'Number of genomes to be downloaded: {str(second_value)}')
    sys.exit(0)

def download_genomes_taxon(genus_species_command,db):
    max_retries = 15

    for attempt in range(max_retries):
        try:
            run_datasets = subprocess.run([f"datasets download genome taxon {genus_species_command} --filename multiple_datasets.zip --assembly-source {db}"], shell=True, capture_output=True, text=True)
            print(run_datasets)
            if run_datasets.returncode == 0:
                logging.info("ncbi-datasets command executed successfully")
                break
            else:
                logging.error(f"Command failed with return code {run_datasets.returncode}")
                if attempt < max_retries - 1:
                    logging.warning("Retrying ncbi-datasets...")
                    time.sleep(1)
                else:
                    if attempt == max_retries:
                        logging.error('Ncbi-datasets not working. Please try again.')
                        sys.exit(1)
        except Exception as e:
            logging.error(f"An error occurred: {str(e)}")

    zip_file = f'multiple_datasets.zip'
    file_extension = '.fna'

    if not os.path.exists(extractedFiles_directory):
        os.mkdir(extractedFiles_directory)

    extract_files_from_zip(zip_file, file_extension, extractedFiles_directory)
    shutil.rmtree(f'{extractedFiles_directory}/ncbi_dataset/')

def download_genomes_accession(accessions_to_download):
    counter = 0
    grouping_factor = 500
    list_of_accessions = [accessions_to_download[i:i+grouping_factor] for i in range(0, len(accessions_to_download), grouping_factor)]
    for lst in list_of_accessions:

        datasets_argument = ' '.join(lst)

        max_retries = 15

        for attempt in range(max_retries):
            try:
                run_datasets = subprocess.run([f"datasets download genome accession {datasets_argument} --filename multiple_datasets_{counter}.zip"], shell=True, capture_output=True, text=True)
                print(run_datasets)
                if run_datasets.returncode == 0:
                    logging.info("ncbi-datasets command executed successfully")
                    break
                else:
                    logging.error(f"Command failed with return code {run_datasets.returncode}")
                    if attempt < max_retries - 1:
                        logging.warning("Retrying ncbi-datasets...")
                        time.sleep(1)
                    else:
                        if attempt == max_retries:
                            logging.error('Ncbi-datasets not working. Please try again.')
                            sys.exit(1)
            except Exception as e:
                logging.error(f"An error occurred: {str(e)}")
        zip_file = f'multiple_datasets_{counter}.zip'
        file_extension = '.fna'

        if not os.path.exists(extractedFiles_directory):
            os.mkdir(extractedFiles_directory)

        extract_files_from_zip(zip_file, file_extension, extractedFiles_directory)
        shutil.rmtree(f'{extractedFiles_directory}/ncbi_dataset/')

        counter += 1

def run_mlst(genome_dict, EXT, sequence_type):
    mlst_command = subprocess.run([f'mlst --nopath {extractedFiles_directory}*.fna > mlst_report.tsv'], shell=True)
    with open('mlst_report.tsv', 'r') as report:
        for line in report:
            line = line.rstrip()
            line_list = line.split()
            genome_name = line_list[0].split(EXT)[0]
            if line_list[2] == sequence_type:
                genome_dict[genome_name] = 'CASE'
            elif line_list[2] == '-':
                logging.info(f'No sequence type scheme available for {genome_name}')
            else:
                genome_dict[genome_name] = 'CONTROL'
    if len(genome_dict) == 0:
        logging.error('An MLST scheme does not exist for this species. Exiting pipeline...')
        sys.exit(1)

def write_case_control_file(TIMESTR, genome_dict):
    output_case_control_file = 'case_control_file_' + TIMESTR + '.txt'
    with open(output_case_control_file, 'w') as out:
        for k,v in genome_dict.items():
            v = v.lower()
            k = k + EXT
            out.write(f'{k},{v}\n')

def create_counter_dictionary(genome_dict):
    value_counter = Counter(genome_dict.values())
    num_controls = value_counter['CONTROL']
    num_cases = value_counter['CASE']

    if (num_controls == 0) or (num_cases == 0):
        logging.error('Either there are no cases or no controls. Exiting pipeline...')
        sys.exit(1)
    return num_controls, num_cases

def generate_smaller_datasets(genome_dict, cases_count, controls_count):
    case_genomes = [genome for genome, status in genome_dict.items() if status == 'CASE']
    control_genomes = [genome for genome, status in genome_dict.items() if status == 'CONTROL']

    genome_lists = []

    while case_genomes or control_genomes:
        selected_cases = case_genomes[:cases_count]
        selected_controls = control_genomes[:controls_count]

        selected_genomes = selected_cases + selected_controls
        genome_lists.append(selected_genomes)

        case_genomes = case_genomes[cases_count:]
        control_genomes = control_genomes[controls_count:]

    return genome_lists

def define_sensitivity_and_specificity(sensitivity, specificity):
    sensitivity = sensitivity / 100
    sensitivity = round(sensitivity, 2)
    specificity = specificity / 100
    specificity = round(specificity, 2)

    sensitivity_step1 = num_cases * sensitivity # 1000 * .9 = 900
    sensitivity_step2 = int(num_cases - sensitivity_step1) # 1000 - 900 = 100

    # Find maximum nuber of controls that can be found
    specificity_step1 = num_controls * specificity # 70000 * .9 = 63,0000
    specificity_step2 = int(num_controls - specificity_step1) # 70,000 - 63,000 = 7000

    logging.info('Minimum number of cases: ' + str(int(sensitivity_step1)))
    logging.info('Maximum number of controls allowed: ' + str(specificity_step2))

    global_sensitivity_threshold = sensitivity_step2 # 100
    global_specificity_threshold = specificity_step2 # 7100

    return global_sensitivity_threshold, global_specificity_threshold, specificity_step2, sensitivity_step1

def use_datasets(case_control, EXT, genome_dict):
    accessions_to_download = []
    for line in open(case_control, 'r', encoding='utf-8-sig'):
        line = line.rstrip()
        line_list = line.split(',')
        genome = line_list[0]
        genome = genome.split(EXT)[0]
        classification = line_list[1]
        if classification == 'case':
            genome_dict[genome] = 'CASE'
        else:
            if classification == 'control':
                genome_dict[genome] = 'CONTROL'
        accessions_to_download.append(genome)
    return accessions_to_download

def case_control_file(case_control,EXT,genome_dict):
    for line in open(case_control, 'r', encoding='utf-8-sig'):
        line = line.rstrip()
        line_list = line.split(',')
        genome = line_list[0]
        genome = genome.split(EXT)[0]
        classification = line_list[1]
        if classification == 'case':
            genome_dict[genome] = 'CASE'
        else:
            if classification == 'control':
                genome_dict[genome] = 'CONTROL'

def dataset_composition(sublist, genome_dict):
    num_cases_in_dataset = 0
    num_controls_in_dataset = 0

    for GENOME in sublist:
        if genome_dict[GENOME] == 'CASE':
            num_cases_in_dataset += 1
        else:
            if genome_dict[GENOME] == 'CONTROL':
                num_controls_in_dataset += 1

    logging.info('Cases: ' + str(num_cases_in_dataset))
    logging.info('Controls: ' + str(num_controls_in_dataset))

    return num_cases_in_dataset, num_controls_in_dataset

def write_refs_file(sublist, EXT, extractedFiles_directory):
    with open('refs.txt', 'w') as input_list:
        for genome in sublist:
            genome = genome + EXT
            if not os.path.exists(f'{extractedFiles_directory}{genome}'):
                logging.error(f'{genome} does not exist. Please make sure that all files are named properly. Exiting...')
                sys.exit(1)
            genome_filename = genome
            input_list.write(f'{extractedFiles_directory}{genome_filename}\n')
    input_list.close()

def calculate_local_thresholds(num_cases_in_dataset, global_sensitivity_threshold, global_specificity_threshold):
    local_sensitivity = num_cases_in_dataset - global_sensitivity_threshold
    if local_sensitivity < 0: # In case the lower limit is a negative number.
        local_sensitivity = 0
    logging.info('Local sensitivity: ' + str(local_sensitivity))
    local_specificity = num_cases_in_dataset + global_specificity_threshold
    logging.info('Local specificity: ' + str(local_specificity))
    return local_sensitivity, local_specificity

def generate_unitig_caller_cmd(unitig_counter, thr, fop, master_files_lst, kmer_length):
    if unitig_counter == 0:
        if ARGS.use_simple:
            print('First iteration - Running simple mode')
            unitig_caller_cmd = 'unitig-caller --simple --refs refs.txt --threads {} --unitigs {} --out {}'.format(thr, kmer_list, fop)
            master_files_lst.append('NA')
        else:
            unitig_caller_cmd = 'unitig-caller --call --refs refs.txt --threads {} --kmer {} --out {}'.format(thr,kmer_length,fop)
            master_files_lst.append('NA')
    else:
        print('Running simple mode')
        master_file = master_files_lst[unitig_counter]
        unitig_caller_cmd = 'unitig-caller --simple --refs refs.txt --threads {} --unitigs {} --out {}'.format(thr,master_file,fop)
    logging.info(unitig_caller_cmd)
    return unitig_caller_cmd

def check_num_lines(fop,unitig_caller):
    wc_cmd = 'wc -l {}.pyseer'.format(fop)
    wc_command = subprocess.run(wc_cmd, shell=True, capture_output=True, text=True)
    try:
        number_of_lines = int(wc_command.stdout.split()[0])
    except:
        logging.warning('warning for {}'.format(fop))
        number_of_lines = -1
    fop = fop +'.pyseer'
    # Check if the number of lines is zero
    if unitig_counter == 0:
        if number_of_lines == 0:
            logging.error("File has no lines. Most probably kmer/unitig ran out of memory or no unique kmers found!")
            sys.exit(1)
        else:
            logging.info("Finished running {} batch/es".format(str(unitig_counter)))
    else:
        if number_of_lines == 0:
            logging.warning("Finished running {} batch/es but the output file has no lines. Most probably kmer/unitig ran out of memory or no overlapping kmers!".format(str(unitig_counter)))
        else:
            logging.info("Finished running {} batch/es".format(str(unitig_counter)))
    return fop

def run_unitig_caller(unitig_caller_cmd, unitig_counter):
    unitigcaller = subprocess.run(unitig_caller_cmd, shell=True, capture_output=True)

    if unitigcaller.returncode != 0:
        logging.error("Error occurred while executing unitig command.")
        logging.error(unitigcaller.stderr.decode())
        logging.error('Exiting script. Unitig-caller did not run properly. Please check log file.')
        sys.exit(1)
    else:
        logging.info('Finished unitig-caller')
    unitig_counter += 1
    return unitig_counter

###############################################################################
if ARGS.summary:
    download_genomes_summary(genus_species_command, db)
###############################################################################
# Initialize genome dictionary
genome_dict = {}
if ARGS.species:

    if ARGS.case_accession_list:

        case_list = ARGS.case_accession_list
        accessions_to_download = []
        for line in open(case_list, 'r'):
            line = line.rstrip()
            accessions_to_download.append(line)

        download_genomes_taxon(genus_species_command, db)

        for file in os.listdir(extractedFiles_directory):
            file = file.split(EXT)[0]
            if file in accessions_to_download:
                genome_dict[file] = 'CASE'
            else:
                genome_dict[file] = 'CONTROL'

        value_to_check = 'CASE'

        if value_to_check in genome_dict.values():
            print('Cases were downloaded')
        else:
            download_genomes_accession(accessions_to_download)

        for file in accessions_to_download:
            genome_dict[file] = 'CASE'

    elif ARGS.case_genomes:

        for file in os.listdir(extractedFiles_directory):
            file = file.split(EXT)[0]
            genome_dict[file] = 'CASE'

        download_genomes_taxon(genus_species_command, db)

        for file in os.listdir(extractedFiles_directory):
            file = file.split(EXT)[0]
            if file in genome_dict:
                continue
            else:
                genome_dict[file] = 'CONTROL'
    else:
        if ARGS.sequence_type:
            download_genomes_taxon(genus_species_command, db)
            run_mlst(genome_dict, EXT, sequence_type)
            #os.remove('mlst_report.tsv')

    write_case_control_file(TIMESTR, genome_dict)
################################################################################
if ARGS.case_control_file and not ARGS.use_datasets:
    case_control_file(case_control,EXT,genome_dict)
################################################################################
if ARGS.use_datasets:
    accessions_to_download = use_datasets(case_control, EXT, genome_dict)
    download_genomes_accession(accessions_to_download)
    write_case_control_file(TIMESTR, genome_dict)
################################################################################
print(genome_dict)
num_controls, num_cases = create_counter_dictionary(genome_dict)

if ARGS.number_of_cases:
    cases_count = ARGS.number_of_cases
else:
    cases_count = 400

if ARGS.number_of_controls:
    controls_count = ARGS.number_of_controls
else:
    controls_count = 200
###############################################################################
if ARGS.sensitivity or ARGS.sensitivity == 0:
    sensitivity = ARGS.sensitivity
else:
    sensitivity = 100
if ARGS.specificity or ARGS.specificity == 0:
    specificity = ARGS.specificity
else:
    specificity = 100

global_sensitivity_threshold, global_specificity_threshold, specificity_step2, sensitivity_step1 = define_sensitivity_and_specificity(sensitivity, specificity)

genome_lists = generate_smaller_datasets(genome_dict, cases_count, controls_count)

master_dictionary = {}
removed_kmers = set()
unitig_counter = 0
master_files_lst = []
simple_mode_kmers = []

# When running simple mode, if the k-mer is not found in the iteration, it needs to still get carried over to search future iterations.
if ARGS.use_simple:
    with open(kmer_list, 'r') as opened_file:
        filelines = opened_file.readlines()
        for line in filelines[1:]:
            line = line.rstrip()
            simple_mode_kmers.append(line)

for sublist in genome_lists:
    num_cases_in_dataset, num_controls_in_dataset = dataset_composition(sublist, genome_dict)
    write_refs_file(sublist, EXT, extractedFiles_directory)

    # Run unitig-caller. Unzip the output.
    TIMESTR = time.strftime("%Y%m%d_%H%M%S")
    fop = "unitig_caller_{}_{}".format(str(unitig_counter),TIMESTR)
    fop = os.path.join(temp_dir, fop)
    unitig_caller_cmd = generate_unitig_caller_cmd(unitig_counter, thr, fop, master_files_lst, kmer_length)
    unitig_counter = run_unitig_caller(unitig_caller_cmd, unitig_counter)
    fop = check_num_lines(fop,unitig_counter)

    # Calculate local sensitivity and specificity thresholds. Default range for searching list length.
    local_sensitivity, local_specificity = calculate_local_thresholds(num_cases_in_dataset, global_sensitivity_threshold, global_specificity_threshold)

    # Parse the unitig-caller output
    try:
        for line in open(fop, 'r'):
            case_counter = 0
            control_counter = 0
            line_list = line.rstrip().split()
            kmer = line_list[0]
            genome_list = line_list[2:]
            if kmer in removed_kmers:
                continue

            elif kmer in master_dictionary:
                existing_remaining_case_count = master_dictionary[kmer][0] # Get the existing number of cases this kmer is found in
                existing_control_count = master_dictionary[kmer][1] # Get the existing number of controls this kmer is found in
                recalc_local_sensitivity = num_cases_in_dataset - existing_remaining_case_count
                if recalc_local_sensitivity < 0: # In case the lower limit is a negative number.
                    recalc_local_sensitivity = 0
                remaining_control_count = global_specificity_threshold - existing_control_count
                recalc_local_specificity = num_cases_in_dataset + remaining_control_count
                if (len(genome_list) >= recalc_local_sensitivity) and (len(genome_list) <= recalc_local_specificity): # If the kmer exists in the master dictionary, need to count number of cases and controls and add them to existing counts in the dictionary. Remove the kmer if necessary
                    for genome in genome_list:
                        genome_name = genome.split(':')[0]
                        genome_name = genome_name
                        if genome_dict[genome_name] == 'CONTROL':
                            control_counter += 1
                        else:
                            if genome_dict[genome_name] == 'CASE':
                                case_counter += 1
                    num_of_cases_kmer_not_found_in = num_cases_in_dataset - case_counter
                    updated_case_count = existing_remaining_case_count - num_of_cases_kmer_not_found_in
                    updated_control_count = existing_control_count + control_counter # Update the number of controls
                    if (updated_case_count < 0) or (updated_control_count > global_specificity_threshold): # Check to make sure all of the cases or controls are not 'spent'. If they are spent, remove the kmer from the master dictionary and add to removed kmers set.
                        del master_dictionary[kmer]
                        removed_kmers.add(kmer)

                    else: # If the cases/controls are not spent, update the dictionary with the new counts.
                        master_dictionary[kmer][0] = updated_case_count
                        master_dictionary[kmer][1] = updated_control_count
                        master_dictionary[kmer][2] += case_counter

                else: # If the kmer exists, but doesn't pass the threshold, remove it from dictionary and add to removed kmers set.
                    del master_dictionary[kmer]
                    removed_kmers.add(kmer)

            elif (len(genome_list) >= local_sensitivity) and (len(genome_list) <= local_specificity): # Execute this block of code if the kmer has not been seen yet, but passes the initial threshold.
                for genome in genome_list:
                    genome_name = genome.split(':')[0]
                    genome_name = genome_name
                    if genome_dict[genome_name] == 'CONTROL':
                        control_counter += 1
                    else:
                        if genome_dict[genome_name] == 'CASE':
                            case_counter += 1

                if control_counter > specificity_step2: # If controls is greater than 7900, automatically fails. Checking to make sure that the number of controls in not over the allowable limit. Needs to be done when adding a kmer to master dictionary the first time.
                    removed_kmers.add(kmer)
                else:
                    num_of_cases_kmer_not_found_in = num_cases_in_dataset - case_counter
                    updated_case_count = global_sensitivity_threshold - num_of_cases_kmer_not_found_in # 100 - number of cases for this kmer
                    if (updated_case_count < 0):
                        removed_kmers.add(kmer)
                    else: # If not spent, add to the master dictionary
                        subdict_value = [updated_case_count, control_counter, case_counter]#check this
                        master_dictionary[kmer] = subdict_value

            else:
                removed_kmers.add(kmer)
    except:
        logging.warning("Error occurred while parsing unitig output.".format(unitig_counter))

    master_file_iter = "{}/unitig_caller_master_{}_{}.txt".format(temp_dir,str(unitig_counter),TIMESTR)
    master_files_lst.append(master_file_iter)
    mast_iter = open(master_file_iter, 'w')
    mast_iter.write('kmer\n')
    mast_iter.writelines('{}\n'.format(k) for k in master_dictionary)
    # When running simple mode, if the k-mer is not found in the iteration, it needs to still get carried over to search future iterations.
    if len(simple_mode_kmers) != 0:
        for kmer in simple_mode_kmers:
            if kmer in removed_kmers:
                continue
            elif kmer in master_dictionary:
                continue
            else:
                mast_iter.write(f'{kmer}\n')
    mast_iter.close()


kmers_to_be_removed = []
for k in master_dictionary:
    if master_dictionary[k][-1] < sensitivity_step1:
        kmers_to_be_removed.append(k)

for k in kmers_to_be_removed:
    del master_dictionary[k]
    removed_kmers.add(k)

with open('UniqueKmers.txt', 'w') as out:
    for k in master_dictionary:
        out.write(k + '\n')

logging.info('Number of removed kmers: ' + str(len(removed_kmers)))
rem_op = open(f'{temp_dir}/removed_kmers.txt','w')
rem_op.writelines(f"{k}\n" for k in removed_kmers)
mast_op = open('Unique_Kmers_Report.txt','w')
mast_op.write('Kmer\tNumber of Controls K-mer Found in\tNumber of Cases K-mer Found In\n')
for k,v in master_dictionary.items():
    kmer = k
    found_in_controls = v[1]
    found_in_cases = v[2]
    mast_op.write(f'{kmer}\t{found_in_controls}\t{found_in_cases}\n')
logging.info('Number of kmers in master_dict: ' + str(len(master_dictionary)))
# os.remove('refs.txt')
