import os
import sys
import csv
import tempfile
import shutil
import pytest

#pytest fixture for tests
# run all tests with "pytest tests/""

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def pytest_configure(config):
    config.addinivalue_line("markers", "slow: marks tests as slow-running")

@pytest.fixture
def temp_dir():
    """
    Creates a temporary direction for the test.
    the 'yeild' returns the path to the test, then cleanup after.
    """
    temp_path = tempfile.mkdtemp()
    yield temp_path #this is the test location
    shutil.rmtree(temp_path) #this cleans up the test enviormoment after the test

@pytest.fixture
def crispresso_output_dir(temp_dir):
    """
    Creates the directory structure that the script expects"
        temp_dir/
            CRISPResso_on_test/
    Returns both paths so tests can put files in the right place.
    """
    crispresso_subdir = os.path.join(temp_dir, "CRISPResso_on_test")
    os.makedirs(crispresso_subdir)
    return temp_dir, crispresso_subdir

class CRISPRessoFileFactory:
    """
    Factory class for creating mock CRISPResso output files.
    """
    @staticmethod
    def create_quant_window_file(crispresso_dir, base_percentages):
        """
        Creates an equivilent to Quantification_window_nucleotide_percentage_table.txt
        
                Base    1       2       3       4       5       ...
        A       0.95    0.02    0.50    0.02    0.95    ...
        C       0.02    0.95    0.50    0.95    0.02    ...
        G       0.02    0.02    0.00    0.02    0.02    ...
        T       0.01    0.01    0.00    0.01    0.01    ...

        returns:
            str: Path to created file

        """
        quant_file = os.path.join(
            crispresso_dir,
            "Quantification_window_nucleotide_percentage_table.txt"
        )
        num_positions = len(base_percentages['A'])
        header = "Base\t" + "\t".join([str(i + 1) for i in range(num_positions)])

        with open(quant_file, 'w', encoding='utf-8-sig') as f:
            f.write(header + "\n")
            for base in ['A', 'C', 'G', 'T']:
                row = base + "\t" + "\t".join([str(p) for p in base_percentages[base]])
                f.write(row + "\n")
        
        return quant_file

    @staticmethod
    def create_alleles_frequency_table(crispresso_dir, rows, filename=None):
        """
        Creates a mock Alleles_frequency_table.txt
        """
        if filename is None:
            filename = "Alleles_frequency_table.txt"

        alleles_file = os.path.join(crispresso_dir, filename)
        header = "Aligned_Sequence\tReference_Sequence\tUnedited\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads"

        total_reads = sum(row[1] for row in rows)

        with open(alleles_file, 'w', encoding='utf-8-sig') as f:
            f.write(header + "\n")
            for seq, n_reads in rows:
                pct = (n_reads / total_reads * 100) if total_reads > 0 else 0
                f.write(f"{seq}\t{seq}\tFalse\t0\t0\t0\t{n_reads}\t{pct:.4f}\n")
        
        return alleles_file

@pytest.fixture
def file_factory():
    """Provides access to the file factory class."""
    return CRISPRessoFileFactory