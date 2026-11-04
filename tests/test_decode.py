import os
import subprocess
import pytest
import shutil

# Get the directory this test script lives in
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Define the relative path to the test_output directory from the test file

TEST_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "./decode/test_output")
PRE_CREATED_FILES_DIR = os.path.join(SCRIPT_DIR, "./decode")

# Define file names to compare
FILE_NAME_1 = 'example_8-2_barcodes_decoded.txt'
FILE_NAME_2 = 'example_8-2_barcodes_prefixes_TCTACTCTCCATACG_CACTTGGATC_decoded.txt'

barcode_file_path = os.path.join(SCRIPT_DIR, "../barcodes/barcodes8-2.txt")
barcode_fq_file_path = os.path.join(SCRIPT_DIR, "../examples/example_8-2_barcodes.fq")
prefixes_file_path = os.path.join(SCRIPT_DIR, "../examples/example_8-2_barcodes_prefixes_TCTACTCTCCATACG_CACTTGGATC.fq")

# Fixture to create and clean up the test directory
@pytest.fixture(scope="module")
def setup_and_teardown():
    # Ensure the test_output directory does not already exist
    if os.path.exists(TEST_OUTPUT_DIR):
        shutil.rmtree(TEST_OUTPUT_DIR)

    # Create the test_output directory
    os.makedirs(TEST_OUTPUT_DIR)

    yield  # Test execution will happen after this point

    # Cleanup: Remove the test_output directory after tests run
    if os.path.exists(TEST_OUTPUT_DIR):
        shutil.rmtree(TEST_OUTPUT_DIR)


# Test function to run the freebarcodes commands
def test_freebarcodes_decode(setup_and_teardown):
    # First command
    command1 = [
        'freebarcodes', 'decode', barcode_file_path,
        barcode_fq_file_path,
        '--output-dir=' + TEST_OUTPUT_DIR
    ]
    
    # Second command
    command2 = [
        'freebarcodes', 'decode', barcode_file_path,
        prefixes_file_path,
        '--prefixes=TCTACTCTCCATACG,CACTTGGATC', '--max-prefix-err=3,2',
        '--output-dir=' + TEST_OUTPUT_DIR
    ]
    # Print the commands to ensure they are correct
    print(f"apple Command 1: {' '.join(command1)}")
    print(f"apple Command 2: {' '.join(command2)}")

    # Run both commands using subprocess
    try:
        # Running the first command
        result1 = subprocess.run(command1, check=True, capture_output=True, text=True)
        
        print(f"Command 1 output:\n{result1.stdout}")
        print(f"Command 1 error:\n{result1.stderr}")

        # Running the second command
        result2 = subprocess.run(command2, check=True, capture_output=True, text=True)
        
        print(f"Command 2 output:\n{result2.stdout}")
        print(f"Command 2 error:\n{result2.stderr}")

    except subprocess.CalledProcessError as e:
        pytest.fail(f"Command failed with error: {e.stderr}")

    # Define the paths for both files in the test output directory and pre-created files directory
    output_file1 = os.path.join(TEST_OUTPUT_DIR, FILE_NAME_1)
    output_file2 = os.path.join(TEST_OUTPUT_DIR, FILE_NAME_2)

    pre_created_file1 = os.path.join(PRE_CREATED_FILES_DIR, FILE_NAME_1)
    pre_created_file2 = os.path.join(PRE_CREATED_FILES_DIR, FILE_NAME_2)

    # Assert that the files have identical content
    assert compare_files(output_file1, pre_created_file1), f"Content of {FILE_NAME_1} does not match"
    assert compare_files(output_file2, pre_created_file2), f"Content of {FILE_NAME_2} does not match"

# why
def compare_files(file1, file2):
    """
    Compare the content of two files and return True if they are the same, False otherwise.
    """
    try:
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            return f1.read() == f2.read()
    except Exception as e:
        pytest.fail(f"Error comparing files {file1} and {file2}: {e}")
        return False