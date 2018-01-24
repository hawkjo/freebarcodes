# FreeBarcodes

A package for the generation and decoding of FREE divergence error-correcting DNA barcodes, as described in the manuscript:

#### Error-correcting DNA barcodes for next-generation sequencing
*John A. Hawkins, Stephen K. Jones Jr, Ilya J. Finkelstein, and William H. Press*

*(Submitted)*

### Installation

The following instructions should work across platforms, except that installing virtualenv with apt-get is Ubuntu specific. For other platforms, install virtualenv appropriately if desired.

First, clone the repository to a local directory:

```
git clone https://github.com/hawkjo/freebarcodes.git
```

Optionally, you can install into a virtual environment (recommended):

```
sudo apt-get install -y virtualenv
cd freebarcodes
virtualenv envfreebarcodes
. envfreebarcodes/bin/activate
```

Now install required packages and freebarcodes using pip:

```
pip install numpy==1.13.3 && pip install -r requirements.txt && python setup.py install
```

### Usage

```
Usage:
  freebarcodes.py decode       <barcode_files> <fastq_files> [--output-dir=<output_dir>] [--prefixes=<prefixes>] [-v | -vv | -vvv]
  freebarcodes.py generate     <barcode_length> <num_errors> [--output-dir=<output_dir>] [-v | -vv | -vvv]
  freebarcodes.py prune        <raw_barcodes_file> <num_errors> [--output-dir=<output_dir>] [-v | -vv | -vvv]
  freebarcodes.py concatenate  <barcode_files> [--output-dir=<output_dir>] [-v | -vv | -vvv]

Options:
  -h --help     Show this screen.
  --version     Show version.

Commands:
  decode        Decode barcodes in fastq files. Separate file names with commas.
  generate      Generate new sets of FREE barcodes of given length and number correctable errors
  prune         Prune a previous set of barcodes to a subset of valid FREE barcodes
  concatenate   Concatenate and filter multiple sets of barcodes
```

#### Pre-generated barcode lists

Most users will be able to use one of the pre-generated lists of barcodes included in  `freebarcodes/barcodes/` for their experiments. Barcodes are stored in lists according to barcode length and number of errors corrected. For example, barcodes of length 15 which correct up to 2 errors are stored in the file `barcodes15-2.txt`. Each line of the file contains a unique barcode.

#### Decode

Decoding barcodes is as easy as telling `freebarcodes.py` which barcode file(s) to use and which fastq file(s) should be decoded. For concatenated barcodes, list all files of barcodes in order, separated by commas and no spaces. Multiple fastq files with the same barcode list can similarly be comma-delimited. For example,
```
freebarcodes.py  decode  barcodes15-2.txt,barcodes15-2.txt  seqdata1.fastq,seqdata2.fastq
```
will decode the barcodes in `seqdata1.fastq` and `seqdata2.fastq`, where the decoding process looks for two 15 bp, 2-error correcting barcodes at the beginning of each read. 

The `--prefixes` option allows the user to look for barcodes after some known starting sequence, for example, after a primer sequence. In the event of multiple primer sequences being possible, comma-delimiting is also allowed and the observed prefix sequence is also output. Only one prefix sequence is allowed per read, however. Note that in some cases the same could be accomplished by using concatenated barcodes, but this method is more general, preferred when the possible prefixes are few (for example, 1 or 2), and possibly long and/or of different lengths.

The results are output in files names after the fastq files. For example, `seqdata1.fastq` results would be stored in `seqdata1_decoded.txt`. Results are stored in the following tab-delimited format, where brackets indicate output from optional arguments:
```
fastq_read_name 	[prefix] 	barcode1 	[barcode2] 	[...] 	sequence
```

#### Generate

Generating new barcodes can be done via the command line using the same algorithm and filters as used for the pre-generated lists by simply providing the desired barcode length and number of errors to correct. The synthesis- and sequencing-friendly filters used are

* Balanced GC content (40-60%)
* No homopolymer triples (e.g., TTT)
* No triplet self-complementarity 
* No GGC (Illumina error motif)

To generate lists using user-specific filtering criteria, a user comfortable in python can create their own iterator of acceptable barcodes and pass it to a `freebarcodes.generate.FreeDivBarcodeGenerator` object to generate a list of valid FREE barcodes. Alternatively, and more simply, one can use the pruning method described below.

#### Prune

To prune a set of pre-existing barcodes to a subset of valid FREE divergence barcodes, make a file with one barcode per line and pass it `freebarcodes.py prune`. The user must also specify a number of errors to correct.

#### Concatenate

For experiments where large numbers of barcodes are needed, one can concatenate barcodes and decode them with the decode command above. Any set of barcodes can be concatenated with any other set of barcodes, but in general one will want to filter out certain barcodes that do not go well together, for example, those prone to internal hairpin structure. To concatenate any set of barcodes with any other set or sets of barcodes, filtering according to the same sequencing and synthesis filters as used for barcode generation, use the `freebarcodes.py concatenate` command with the lists of barcodes to concatenate, in order, and separated by commas and no spaces. Note that GC content is not checked, however, as it is already presumed balanced as desired.


