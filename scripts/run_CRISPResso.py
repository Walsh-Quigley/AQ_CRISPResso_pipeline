import os
import glob
import subprocess

def run_CRISPResso(guide_seq, amplicon_seq, orientation, editor, directory_path):
    print(f"Running CRISPResso for editor: {editor} in directory: {directory_path}")
    # Look for fastq or fastq.gz files in the directory
    fastq_files = glob.glob(os.path.join(directory_path, "*.fastq")) + \
                  glob.glob(os.path.join(directory_path, "*.fastq.gz"))
    
    if not fastq_files:
        print(f"No FASTQ files found in {directory_path}. Skipping...")
        return
    
    elif len(fastq_files) == 1:
            fastq_cmd_section = ['--fastq_r1', fastq_files[0]]
    elif len(fastq_files) == 2:
            fastq_cmd_section = ['--fastq_r1', fastq_files[0], '--fastq_r2', fastq_files[1]]
    else:
        print(f"More than two FASTQ files found in {directory_path}. Please check the files. Skipping...")
        return
    
    
    cmd = [
        'CRISPResso',
        *fastq_cmd_section,
        '--amplicon_seq', amplicon_seq,
        '--guide_seq', guide_seq,
        '--output_folder', directory_path,
        '--plot_window_size', str((len(guide_seq)+1)//2),
        '--quantification_window_center', str(-(len(guide_seq) +1)//2),
        '--quantification_window_size', str((len(guide_seq)+1)//2),
    ]

    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print("\n\033[1;31mCRISPResso failed\033[0m")
        print(f"Error: {result.stderr}")




