import os
import glob
import subprocess
import logging

def run_CRISPResso(guide_seq, amplicon_seq, orientation, editor, directory_path):
    """
    Run CRISPResso on FASTQ files in the specified directory.

    Returns:
        tuple: (success: bool, error_type: str or None)
            error_type can be: 'empty_fastq', 'spaces_in_path', 'no_fastq', 'too_many_fastq', 'other', or None if success
    """
    print(f"Running CRISPResso for editor: {editor} in directory: {directory_path}")
    # Look for fastq or fastq.gz files in the directory
    original_dir = os.getcwd()
    os.chdir(directory_path)

    try:
        # Look for fastq or fastq.gz files using relative paths
        fastq_files = glob.glob("*.fastq") + glob.glob("*.fastq.gz")
        if not fastq_files:
            print(f"No FASTQ files found in {directory_path}. Skipping...")
            return (False, 'no_fastq')

        elif len(fastq_files) == 1:
            fastq_cmd_section = ['--fastq_r1', fastq_files[0]]
        elif len(fastq_files) == 2:
            fastq_cmd_section = ['--fastq_r1', fastq_files[0], '--fastq_r2', fastq_files[1]]
        else:
            print(f"More than two FASTQ files found in {directory_path}. Please check the files. Skipping...")
            return (False, 'too_many_fastq')
        
        
        cmd = [
            'CRISPResso',
            *fastq_cmd_section,
            '--amplicon_seq', amplicon_seq,
            '--guide_seq', guide_seq,
            '--output_folder', '.',
            '--plot_window_size', str((len(guide_seq)+1)//2),
            '--quantification_window_center', str(-len(guide_seq) //2),
            '--quantification_window_size', str((len(guide_seq)+1)//2),
        ]

        logging.info(cmd)

        # Run CRISPResso and capture output while still displaying it in real-time
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        # Collect output while printing it
        output_lines = []
        for line in process.stdout:
            print(line, end='')  # Print to terminal in real-time
            output_lines.append(line)

        process.wait()
        full_output = ''.join(output_lines)

        # Check for 0 reads error
        if "0 reads" in full_output:
            # Check if the FASTQ files are actually empty
            total_size = sum(os.path.getsize(f) for f in fastq_files)
            if total_size < 1000:  # Less than 1KB means likely empty
                error_msg = f"ERROR: Empty FASTQ file for {directory_path}"
                logging.error(error_msg)
                print(f"\n\033[1;31m{error_msg}\033[0m")
                return (False, 'empty_fastq')
            # If files have data but 0 reads aligned, fall through to normal error handling

        if process.returncode != 0:
            # Check if failure is likely due to spaces in file path
            if " " in directory_path or any(" " in f for f in fastq_files):
                error_msg = f"ERROR: CRISPResso failed - likely due to spaces in file path: {directory_path}"
                logging.error(error_msg)
                logging.error(f"CRISPResso output: {full_output[-500:] if len(full_output) > 500 else full_output}")
                print(f"\n\033[1;31mCRISPResso failed\033[0m")
                return (False, 'spaces_in_path')
            else:
                error_msg = f"CRISPResso failed for {directory_path}"
                logging.error(error_msg)
                logging.error(f"CRISPResso output: {full_output[-500:] if len(full_output) > 500 else full_output}")
                print(f"\n\033[1;31mCRISPResso failed\033[0m")
                return (False, 'other')

        return (True, None)
    finally:
        # Always return to the original directory
        os.chdir(original_dir)




