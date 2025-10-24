import os
import glob
import subprocess

def run_CRISPResso(guide_seq, amplicon_seq, orientation, editor, directory_path):

    # Look for fastq or fastq.gz files in the directory
    fastq_files = glob.glob(os.path.join(directory_path, "*.fastq")) + \
                  glob.glob(os.path.join(directory_path, "*.fastq.gz"))
    
    if not fastq_files:
        print(f"No FASTQ files found in {directory_path}. Skipping...")
        return
    
    elif len(fastq_files) == 1:
        fastq_file = fastq_files[0]

        # First attempt with quantification window
        cmd = [
            'CRISPResso',
            '--fastq_r1', fastq_file,
            '--amplicon_seq', amplicon_seq,
            '--guide_seq', guide_seq,
            '--quantification_window_size', '10',
            '--quantification_window_center', '-10',
            '--base_editor_output',
            '--output_folder', directory_path
        ]

        print(f"Running CRISPResso on {fastq_file} with quantification window...")
        result = subprocess.run(cmd)

        # Retry without window if it fails
        if result.returncode != 0:
            print("\nCRISPResso failed with quantification window. Retrying without window...")
            retry_cmd = [
                'CRISPResso',
                '--fastq_r1', fastq_file,
                '--amplicon_seq', amplicon_seq,
                '--guide_seq', guide_seq,
                '--base_editor_output',
                '--output_folder', directory_path
            ]
            subprocess.run(retry_cmd)
