from pathlib import Path
from config import AmpliconConfig
from glob import glob
import subprocess
import logging


#stage 1 -> finds FASTQs, matches each to an amplicon config, runs CRISPResso

def by_name_length(config) -> int:
    """Sort key - return the length of the amplicon's name"""
    return len(config.name)

def pair_fastq_files(fastq_files: list[str]) -> tuple[str, str]:
    """Identifies R1 and R2 from a list of exactly two fastq file paths.
    Args:
        fastq_files: a list of exactly two fastq file paths
    Returns:
        tuple[str, str]: a tuple of (R1 path, R2 path)
    Raises:
        ValueError: if R1 and R2 cannot be unambiguously identified
    """
    r1_files = [f for f in fastq_files if "_R1" in Path(f).name.upper()]
    r2_files = [f for f in fastq_files if "_R2" in Path(f).name.upper()]
    if len(r1_files) == 1 and len(r2_files) == 1:
        read1 = r1_files[0]
        read2 = r2_files[0]
    else:
        raise ValueError(f"could not unambiguously identify R1/R2 in {fastq_files}")
    return read1, read2

def build_window_args(amplicon_list_row: AmpliconConfig) -> list[str]:
    """checks editor and creates the correct window arguments for the CRISPResso command, 
        based on the editor.
    Args:
        amplicon_list_row: the AmpliconConfig object used for this CRISPResso call
    Returns:
        list[str]: a list of arguments to be passed to the CRISPResso call
    """
    proto_len = len(amplicon_list_row.protospacer)
    if amplicon_list_row.editor == "NUCLEASE":
        return [
            '--quantification_window_center', '-3',
            '--quantification_window_size', '15',
            '--plot_window_size', '15',
            '--default_min_aln_score', str(amplicon_list_row.min_alignment_score),
        ]
    return [
            '--plot_window_size', str((proto_len + 1) // 2),
            '--quantification_window_center', str(-proto_len // 2),
            '--quantification_window_size', str((proto_len + 1) // 2),
        ]

def identify_amplicon(directory_name: str, amplicon_configs: list[AmpliconConfig]) -> AmpliconConfig:
    """matches the correct amplicon to the given sample
    Args:
        directory_name: the name of the sample directory
        amplicon_configs: list of all AmpliconConfig objects from amplicon_list.csv
    Returns:
        AmpliconConfig: the amplicon that matches the sample directory
    Raises:
        ValueError: no amplicon config object matches the sample directory
    """
    
    matched_name = None

    clean_name = directory_name.split(".")[0]

    directory_upper = clean_name.upper()

    for config in sorted(amplicon_configs, key=by_name_length, reverse=True):        
        if config.name.upper() in directory_upper:            
            matched_name = config
            break
    
    if not matched_name:
        error_msg = f"No valid amplicon match found for directory: {directory_name}"
        raise ValueError(error_msg)

    logging.info(f"Matched {directory_name} to amplicon {matched_name.name}")
    return matched_name

def run_crispresso(amplicon_list_row: AmpliconConfig, sample_dir: Path) -> None:
    """Runs the CRISPResso command line function using information from the matched amplicon
        config file
    Args:
        amplicon_list_row: the AmpliconConfig object that was associated with the sample
        sample_dir: the directory path of the current sample analysis is being run on
    Returns:
        None: the purpose of the function is to run the CRISPResso command, no return value
    Raises:
        FileNotFoundError: no fastq files found in the sample directory
        ValueError: unable to distinguish read 1 and read 2 in paired end reads
        ValueError: more than 2 fastq files found in the sample directory
    """
    fastq_files = sorted(glob(str(sample_dir / "*.fastq.gz")) + glob(str(sample_dir / "*.fastq")))

    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files found in {sample_dir}")

    elif len(fastq_files) == 1:
        fastq_cmd_section = ['--fastq_r1', fastq_files[0]]
    elif len(fastq_files) == 2:
        read1, read2 = pair_fastq_files(fastq_files)
        fastq_cmd_section = ['--fastq_r1', read1, '--fastq_r2', read2]
    else:
        raise ValueError(f"More than two FASTQ files found in {sample_dir}")

    ####Static Args the crispresso command need regardless of editor
    common_args = [
        'CRISPResso',
        *fastq_cmd_section,
        '--amplicon_seq', amplicon_list_row.amplicon, #amplicon sequence from amplicon config object
        '--guide_seq', amplicon_list_row.protospacer, #protospacer sequence from amplicon config object
        '--output_folder', str(sample_dir), #output folder for the crispresso run
    ]

    #helper funtion that populates the remaining window args based on editor
    window_args = build_window_args(amplicon_list_row)

    cmd = common_args + window_args
    subprocess.run(cmd, check=True)

