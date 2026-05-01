from pathlib import Path
from config import AmpliconConfig
from glob import glob
import subprocess


#stage 1 -> finds FASTQs, matches each to an amplicon config, runs CRISPResso

def by_name_length(config):
    return len(config.name)


def identify_amplicon(directory_name: str, amplicon_configs: list[AmpliconConfig]) -> AmpliconConfig:
    
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

    return matched_name

def run_crispresso(amplicon_list_row: AmpliconConfig, sample_dir: Path) -> None:
    fastq_files = glob(str(sample_dir / "*.fastq.gz")) + glob(str(sample_dir / "*fastq"))

    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files found in {sample_dir}")

    elif len(fastq_files) == 1:
        fastq_cmd_section = ['--fastq_r1', fastq_files[0]]
    elif len(fastq_files) == 2:
        fastq_cmd_section = ['--fastq_r1', fastq_files[0], '--fastq_r2', fastq_files[1]]
    else:
        raise ValueError(f"More than two FASTQ files found in {sample_dir}")

    cmd = [
        'CRISPResso',
        *fastq_cmd_section,
        '--amplicon_seq', amplicon_list_row.amplicon,
        '--guide_seq', amplicon_list_row.protospacer,
        '--output_folder', str(sample_dir),
        '--plot_window_size', str((len(amplicon_list_row.protospacer) + 1) // 2),
        '--quantification_window_center', str(-len(amplicon_list_row.protospacer) //2),
        '--quantification_window_size', str((len(amplicon_list_row.protospacer)+1)//2),
    ]    

    subprocess.run(cmd, check=True)

