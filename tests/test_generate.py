import os
import subprocess
import pytest
import shutil

# Get the directory this test script lives in
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Define the relative path to the test_output directory from the test file

TEST_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "./generate/test_output")
PRE_CREATED_FILES_DIR = os.path.join(SCRIPT_DIR, "./generate")

# Define file names to compare
FILE_NAME = "barcodes10-2.txt"

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
def test_freebarcodes_generate(setup_and_teardown):
    # First command
    command = [
        'freebarcodes', 'generate', "10", "2",
        '--output-dir=' + TEST_OUTPUT_DIR
    ]

    # Print the commands to ensure they are correct
    print(f"apple Command 1: {' '.join(command)}")

    # Run both commands using subprocess
    try:
        # Running the first command
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        
        print(f"Command 1 output:\n{result.stdout}")
        print(f"Command 1 error:\n{result.stderr}")

    except subprocess.CalledProcessError as e:
        pytest.fail(f"Command failed with error: {e.stderr}")

    # Define the paths for both files in the test output directory and pre-created files directory
    output_file = os.path.join(TEST_OUTPUT_DIR, FILE_NAME)

    pre_created_file = os.path.join(PRE_CREATED_FILES_DIR, FILE_NAME)

    # Assert that the files have identical content
    assert compare_files(output_file, pre_created_file), f"Content of {FILE_NAME} does not match"

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