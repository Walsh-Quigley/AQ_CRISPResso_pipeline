"""
End to end test for the complete CRISPResso pipeline

Tests CRISPResso_Loop -> Quantificiation_loop -> output validation.

Run with "pytest tests/test_end_to_end.py -v -s"
The -s flag shows the print output
"""
import os
import sys
import csv
import gzip
import shutil
import tempfile
import pytest
import pandas as pd
from unittest.mock import patch

# Add project root to path so we can import the main scripts
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def apply_edits_to_sequece(sequence, edits):
    """
    Apply base substitutions to a sequence (0 indexed)
    
    Parameters:
        sequence: The DNA sequence string
        edits: Dict of {position: new_base} where position is 0-indexed

    Returns:
        Modified sequence with edits applied

    Example:
        apply_edits_to_sequence("ACGT", {1:'T'}) 
        'ATGT
    """
    seq_list = list(sequence)
    for pos, new_base in edits.items():
        seq_list[pos] = new_base
    return ''.join(seq_list)

def write_fastq_gz(filepath, reads_data):
    """
    Write reads to a gzipped FASTQ file
    
    Args:
        filepath: Path to ouput .fastq.gz file
        reads_data: List of touples: [(sequence, count), ...]
                        Each sequence will be written 'count' times
    
    The fucntion generates simple FASTQ entries with:
    - Unique read IDs
    - The providede sequence
    - High quality scores so we dont filter out any
    """
    read_num = 0
    with gzip.open(filepath, 'wt') as f:
        for sequence, count in reads_data:
            for _ in range(count):
                #FASTQ format: 4 lines per read
                f.write(f"@read_{read_num}\n")      # Header
                f.write(f"{sequence}\n")            # Sequence
                f.write("+\n")                      # Separator
                f.write("I" * len(sequence) + "\n") # Quality scores
                read_num += 1

def find_guide_in_amplicon(amplicon, guide, orientation):
    """
    Find where the guide sequece appears in the amplicon 

    Args:
        amplicon: Full amplicon sequence
        guide: Guide/protospacer sequence (20bp typically)
        orientation: 'F' for forward, 'R' for reverse
    
    Returns:
        Start position (0-indexed) of guide in amplicon

    for reverse orientation, we need the reverse complement
    """
    if orientation == 'F':
        pos = amplicon.upper().find(guide.upper())
        if pos == -1:
            raise ValueError(f"Guide {guide} not found in amplicon")
        return pos
    
    else:
        from scripts.reverse_complement import reverse_complement
        rc_guide = reverse_complement(guide)
        pos = amplicon.upper().find(rc_guide.upper())
        if pos == -1:
            raise ValueError(f"Reverse complement of guide not found in amplicon")
        return pos
    
TEST_AMPLICONS = {
    #case 1, standard ABE editing, PAH1
    "PAH1": {
        "name": "PAH1",
        "guide": "TCACAGTTCGGGGGTATACA",
        "amplicon": "CCTTTTTTTagatggcgctcattgtgcctggcaactggtagctggaggacagtactcacAGTTCGGGGGTATACATGGGCTTGGATCCATGTCTGATGTACTGTGTGCAGCAAGACCTCAATCCTTTGGGTGTATGGGTCG",
        "orientation": "F",
        "editor": "ABE",
        "intended_edit": "5",
        "tolerated_edits": "3,16",
    },

    #case 2, ONE-seq, novel/fictional sequence
    "TEST_ONESEQ": {
        "name": "TEST_ONESEQ",
        "guide": "GGATGGATGGATGGATGGCC",  # 20bp, A's scattered throughout
        "amplicon": "CTGCCTGCCTGCCTGCCTGCCTGCCTGCCTGCCTGCCTGCGGATGGATGGATGGATGGCCCGTCCGTCCGTCCGTCCGTCCGTCCGTCCGTCCGTCCGTC",
        "orientation": "F",
        "editor": "ABE",
        "intended_edit": "ONE-seq",
        "tolerated_edits": "",
    },

    #case 3, het ABE case, novel/fictional
    "TEST_HET": {
        "name": "TEST_HET",
        "guide": "GGGGAGGGGGAGGGGGGGGG",  # A at position 5 (edit target) and position 11
        "amplicon": "TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTC"
                    "GGGGAGGGGGAGGGGGGGGG"
                    "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        "orientation": "F",
        "editor": "ABE",
        "intended_edit": "5",
        "tolerated_edits": "",
    },
}

EXPECTED_RESULTS = {
    # PAH1: Standard ABE case
    # Total: 1000 reads
    "PAH1": {
        "total_reads": 1000,
        "expected_correction_without_bystanders": 30.0,       # 300/1000
        "expected_correction_with_bystanders": 50.0,          # (300+100+100)/1000
        "expected_correction_with_any_AtoG_change": 60.0,     # (300+100+100+100)/1000
        "expected_correction_with_any_change_in_protospacer": 70.0,  # all corrected reads
        "tolerance": 2.0,
    },
    # TEST_ONESEQ: ONE-seq case
    # Total: 1000 reads
    # Mix: 500 no edits, 300 A→G in first 10bp, 200 A→G outside first 10bp
    "TEST_ONESEQ": {
        "total_reads": 1000,
        "expected_AtoG_first_10bp": 30.0,      # 300/1000
        "expected_AtoG_protospacer": 50.0,     # (300+200)/1000
        "tolerance": 2.0,
    },
    
    # TEST_HET: Heterozygous ABE case
    # Total: 1000 reads (500 allele1, 500 allele2)
    "TEST_HET": {
        "total_reads": 1000,
        "expected_correction_without_bystanders": 20.0,  # 200/1000 overall
        "expected_allele1_correction": 40.0,             # 200/500
        "expected_allele2_correction": 0,                # 0/500
        "expected_het_position": 11,                     # 1-indexed position in guide
        "expected_het_alleles": "A/C",
        "tolerance": 3.0,  # Slightly more lenient for het complexity
    },
}

def create_pah1_read_variants():
    """
    Create all sequence variants for PAH1 test case
    
    return a dict with 6 variants to test all 4 correction metrics
    """

    amplicon = TEST_AMPLICONS["PAH1"]["amplicon"].upper()
    guide = TEST_AMPLICONS["PAH1"]["guide"]
    guide_start = find_guide_in_amplicon(amplicon, guide, "F")
    
    # Positions in amplicon (0-indexed)
    pos5 = guide_start + 4    # Intended edit position (A→G)
    pos3 = guide_start + 2    # Tolerated bystander
    pos16 = guide_start + 15  # Tolerated bystander
    pos8 = guide_start + 17    # NON-tolerated A position (for any_AtoG test)
    pos12 = guide_start + 11  # G position for non-A→G change test

    return {
        'uncorrected': amplicon,
        'corrected_only': apply_edits_to_sequece(
            amplicon, {pos5: 'G'}
        ),
        'corrected_plus_bystander3': apply_edits_to_sequece(
            amplicon, {pos5: 'G', pos3: 'G'}
        ),
        'corrected_plus_bystander16': apply_edits_to_sequece(
            amplicon, {pos5: 'G', pos16: 'G'}
        ),
        # Corrected + A→G at non-tolerated position
        'corrected_plus_nontolerated_AtoG': apply_edits_to_sequece(
            amplicon, {pos5: 'G', pos8: 'G'}
        ),
        # Corrected + non-A→G change (G→T)
        'corrected_plus_other_change': apply_edits_to_sequece(
            amplicon, {pos5: 'G', pos12: 'T'}
        ),
    }

def create_oneseq_read_variants():
    """
    Create sequence variants for TEST_ONESEQ case
    
    Guide: GGATGGATGGATGGATGGCC
    A positions: 3, 7, 11, 15 (0-indexed: 2, 6, 10, 14)
    First 10bp: positions 0-9, so A's at indices 2 and 6 are in first 10bp
    """
    amplicon = TEST_AMPLICONS["TEST_ONESEQ"]["amplicon"]
    guide = TEST_AMPLICONS["TEST_ONESEQ"]["guide"]
    guide_start = find_guide_in_amplicon(amplicon, guide, "F")
    
    pos_A_first10 = guide_start + 2 
    pos_A_outside = guide_start + 10
    
    return {
        'no_edits': amplicon,
        
        # A→G in first 10bp of guide
        'AtoG_first10bp': apply_edits_to_sequece(
            amplicon, {pos_A_first10: 'G'}
        ),
        
        # A→G outside first 10bp (but in protospacer)
        'AtoG_outside10bp': apply_edits_to_sequece(
            amplicon, {pos_A_outside: 'G'}
        ),
    }

def create_het_read_variants():
    """
    Create sequence variants for TEST_HET case.
    
    Guide: GGGGAGGGGGAGGGGGGGGG
    - Position 5 (index 4): intended edit A→G
    - Position 11 (index 10): het SNP site (A vs C)
    
    We create two "alleles" that differ at position 11.
    """
    amplicon = TEST_AMPLICONS["TEST_HET"]["amplicon"]
    guide = TEST_AMPLICONS["TEST_HET"]["guide"]
    guide_start = find_guide_in_amplicon(amplicon, guide, "F")
    
    # Positions in amplicon (0-indexed)
    pos5 = guide_start + 4    # Intended edit (A→G)
    pos11 = guide_start + 10  # Het SNP site (A vs C)
    
    # Allele 1: original (has A at position 11)
    allele1_uncorrected = amplicon
    allele1_corrected = apply_edits_to_sequece(amplicon, {pos5: 'G'})
    
    # Allele 2: has C at position 11 (the het SNP)
    allele2_uncorrected = apply_edits_to_sequece(amplicon, {pos11: 'C'})
    allele2_corrected = apply_edits_to_sequece(amplicon, {pos5: 'G', pos11: 'C'})
    
    return {
        'allele1_uncorrected': allele1_uncorrected,
        'allele1_corrected': allele1_corrected,
        'allele2_uncorrected': allele2_uncorrected,
        'allele2_corrected': allele2_corrected,
    }

def create_pah1_samples(fastqs_dir):
    """
    Create PAH1 sample directory with synthetic FASTQ
    """
    sample_dir = os.path.join(fastqs_dir, "TestSample_PAH1")
    os.makedirs(sample_dir)
    
    variants = create_pah1_read_variants()
    
    reads_data = [
        (variants['uncorrected'], 300),
        (variants['corrected_only'], 300),
        (variants['corrected_plus_bystander3'], 100),
        (variants['corrected_plus_bystander16'], 100),
        (variants['corrected_plus_nontolerated_AtoG'], 100),
        (variants['corrected_plus_other_change'], 100),
    ]

    fastq_path = os.path.join(sample_dir, "test_R1.fastq.gz")
    write_fastq_gz(fastq_path, reads_data)

def create_oneseq_sample(fastqs_dir):
    """
    Create TEST_ONESEQ sample directory with sythetic FASTQ.
    """
    sample_dir = os.path.join(fastqs_dir, "TestSample_TEST_ONESEQ")
    os.makedirs(sample_dir)

    variants = create_oneseq_read_variants()

    reads_data = [
        (variants['no_edits'], 500),
        (variants['AtoG_first10bp'], 300),
        (variants['AtoG_outside10bp'], 200),
    ]
    
    fastq_path = os.path.join(sample_dir, "test_R1.fastq.gz")
    write_fastq_gz(fastq_path, reads_data)

def create_het_sample(fastqs_dir):
    """
    Create TEST_HET sample directory with synthetic FASTQ.
    """
    sample_dir = os.path.join(fastqs_dir, "TestSample_TEST_HET")
    os.makedirs(sample_dir)
    
    variants = create_het_read_variants()
    
    reads_data = [
        (variants['allele1_uncorrected'], 300),
        (variants['allele1_corrected'], 200),
        (variants['allele2_uncorrected'], 350),
        (variants['allele2_corrected'], 150),
    ]
    
    fastq_path = os.path.join(sample_dir, "test_R1.fastq.gz")
    write_fastq_gz(fastq_path, reads_data)

def write_test_amplicon_csv(filepath):
    """
    creates a test amplicon_list.csv for our test cases
    """
    header = ["name","protospacer_or_PEG","editor","guide_orientation_relative_to_amplicon","amplicon","note","tolerated_edits","intended_edit"]

    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for name, data in TEST_AMPLICONS.items():
            writer.writerow([
                data["name"],
                data["guide"],
                data["editor"],
                data["orientation"],
                data["amplicon"],
                f"test_{name}",
                data["tolerated_edits"],
                data["intended_edit"],
            ])

def assert_within_tolerance(actual, expected, tolerance, metric_name):
    """
    Assert that actual value is within tolerance of expected value
    """
    diff = abs(actual - expected)
    assert diff <= tolerance, \
        f"{metric_name}: expected {expected}±{tolerance}, got {actual} (diff={diff:.2f})"
    print(f"[TEST]   {metric_name}: {actual:.1f}% (expected {expected}±{tolerance})")




@pytest.fixture
def e2e_test_enviroment():
    """
    Creates a complete isolated test enviroment
    1) temp directory with fatsqs/ subdirectory
    2) test amplicon_list.csv
    3) sample directories with sythetic FASTQ files

    returns enviroment info, then cleans up after test completes.
    """
    #original working direcotry
    original_cwd = os.getcwd()
    #temp direcotry
    temp_dir = tempfile.mkdtemp(prefix="crispresso_e2e_")

    try:
        #create subdirectories
        fastqs_dir = os.path.join(temp_dir, "fastqs")
        os.makedirs(fastqs_dir)

        #log direcotry
        logs_dir = os.path.join(temp_dir, "logs")
        os.makedirs(logs_dir)

        #write test amplicon_list.csv
        amplicon_csv_path = os.path.join(temp_dir, "amplicon_list.csv")
        write_test_amplicon_csv(amplicon_csv_path)

        #create all sampe dircetories with fastqs
        create_pah1_samples(fastqs_dir)
        create_oneseq_sample(fastqs_dir)
        create_het_sample(fastqs_dir)

        #Yield enviroment info to the test
        yield {
            'temp_dir': temp_dir,
            'fastqs_dir': fastqs_dir,
            'original_cwd': original_cwd,
            'amplicon_csv': amplicon_csv_path,
        }

    finally:
        #restore original direcotry and clean up
        os.chdir(original_cwd)
        shutil.rmtree(temp_dir, ignore_errors = True)


@pytest.mark.slow
def test_full_pipeline(e2e_test_enviroment):
    """
    test purpose:
    1) changes to the test directory
    2) runs CRISPResso_Loop to process FASTQs
    3) runs Quantification_loop to generate summaries
    4) verifies output CSV files have expected values
    """
    env = e2e_test_enviroment

    #STEP 1: changing to test directory
    os.chdir(env['temp_dir'])
    print(f"\n[TEST] Working directory: {env['temp_dir']}")
    print(f"[TEST] Amplicon CSV: {env['amplicon_csv']}")
    
    #STEP 2: run the CRISPResso loop
    print("\n[TEST] Running CRISPResso_Loop...")
    from CRISPResso_Loop import main as crispresso_main
    crispresso_main()

    for sample_name in ["TestSample_PAH1", "TestSample_TEST_ONESEQ", "TestSample_TEST_HET"]:
        sample_path = os.path.join(env['fastqs_dir'], sample_name)
        crispresso_dirs = [d for d in os.listdir(sample_path) if d.startswith("CRISPResso_on_")]
        assert len(crispresso_dirs) > 0, f"No CRISPResso output found for {sample_name}"
        print(f"[TEST] CRISPResso output found for {sample_name}")

        #verify alleles frequency table actually exists
        crispresso_path = os.path.join(sample_path, crispresso_dirs[0])
        allele_files = [f for f in os.listdir(crispresso_path) if f.startswith("Alleles_frequency_table")]
        assert len(allele_files) > 0, f"CRISPResso failed for {sample_name} - no Alleles_frequency_table"
        print(f"[TEST] CRISPResso output found for {sample_name}")
    
    #reset working directory back to temp dir before quant loop
    os.chdir(env['temp_dir'])


    #STEP 3: Run Quantification Loop
    print("\n[TEST] Running Quantification_loop...")

    #mock the yes no prompt so we dont get stuck on the prism csv generation
    with patch('scripts.yes_no.yes_no', return_value=False):
        from Quantification_loop import main as quant_main
        quant_main()
    
    #STEP 4: Verify output files exist
    quant_summary = os.path.join(env['fastqs_dir'], "quantification_summary.csv")
    oneseq_summary = os.path.join(env['fastqs_dir'], "quantification_summary_ONE-seq.csv")

    assert os.path.exists(quant_summary), "quantification_summary.csv not created"
    assert os.path.exists(oneseq_summary), "quantification_summary_ONE-seq.csv not created"
    print("[TEST] Output CSV files created")

    #STEP 5: Verify PAH1 results
    df = pd.read_csv(quant_summary)
    pah1_row = df[df['sample'].str.contains('PAH1', case=False)]

    assert len(pah1_row) == 1, f"Expected 1 PAH1 row, found {len(pah1_row)}"
    pah1 = pah1_row.iloc[0]
    expected = EXPECTED_RESULTS["PAH1"]

    # Check correction metrics
    assert_within_tolerance(
        pah1['correction_without_bystanders'],
        expected['expected_correction_without_bystanders'],
        expected['tolerance'],
        "PAH1 correction_without_bystanders"
    )
    assert_within_tolerance(
        pah1['correction_with_bystanders'],
        expected['expected_correction_with_bystanders'],
        expected['tolerance'],
        "PAH1 correction_with_bystanders"
    )
    assert_within_tolerance(
        pah1['correction_with_any_AtoG_change'],
        expected['expected_correction_with_any_AtoG_change'],
        expected['tolerance'],
        "PAH1 correction_with_any_AtoG_change"
    )
    assert_within_tolerance(
        pah1['correction_with_any_change_in_protospacer'],
        expected['expected_correction_with_any_change_in_protospacer'],
        expected['tolerance'],
        "PAH1 correction_with_any_change_in_protospacer"
    )
    print("[TEST] PAH1 results verified")

    #STEP 6: Verify ONE-seq results
    df_oneseq = pd.read_csv(oneseq_summary)
    oneseq_row = df_oneseq[df_oneseq['sample'].str.contains('TEST_ONESEQ', case = False)]

    assert len(oneseq_row) == 1, f"Expected 1 TEST_ONESEQ row, found {len(oneseq_row)}"
    oneseq = oneseq_row.iloc[0]
    expected_oneseq = EXPECTED_RESULTS["TEST_ONESEQ"]

    assert_within_tolerance(
        oneseq['Percent_of_reads_with_A>G_in_first_10bp'],
        expected_oneseq['expected_AtoG_first_10bp'],
        expected_oneseq['tolerance'],
        "ONE-seq A>G in first 10bp"
    )
    assert_within_tolerance(
        oneseq['Percent_of_reads_with_A>G_in_protospacer'],
        expected_oneseq['expected_AtoG_protospacer'],
        expected_oneseq['tolerance'],
        "ONE-seq A>G in protospacer"
    )
    print("[TEST] ONE-seq results verified")

    #STEP 7: Verify HET results
    het_row = df[df['sample'].str.contains('TEST_HET', case=False)]

    assert len(het_row) == 1, f"Expected 1 TEST_HET row, found {len(het_row)}"
    het = het_row.iloc[0]
    expected_het = EXPECTED_RESULTS["TEST_HET"]

    # Check that het was actually detected
    assert pd.notna(het.get('het_position')), "Het position not detected for TEST_HET"
    print(f"[TEST] Het detected at position {het['het_position']}")

    # Check overall correction
    assert_within_tolerance(
        het['correction_without_bystanders'],
        expected_het['expected_correction_without_bystanders'],
        expected_het['tolerance'],
        "TEST_HET overall correction"
    )

    # Check per-allele correction
    assert_within_tolerance(
        het['correction_without_bystanders_allele1'],
        expected_het['expected_allele1_correction'],
        expected_het['tolerance'],
        "TEST_HET allele1 correction"
    )
    assert_within_tolerance(
        het['correction_without_bystanders_allele2'],
        expected_het['expected_allele2_correction'],
        expected_het['tolerance'],
        "TEST_HET allele2 correction"
    )
    print("[TEST] TEST_HET results verified")

    print("\n[TEST] *** ALL TESTS PASSED *** (Yippe!) ")
