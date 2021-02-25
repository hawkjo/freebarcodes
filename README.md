# FreeBarcodes

A package for the generation and decoding of FREE divergence error-correcting DNA barcodes, as described in the manuscript:

### Indel-correcting DNA barcodes for high-throughput sequencing

**John A Hawkins, Stephen K Jones Jr, Ilya J Finkelstein, and William H Press**

*Proc Natl Acad Sci*. June 20, 2018. doi:10.1073/pnas.1802640115.

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
  freebarcodes decode       <barcode_files> <fastq_files> [--output-dir=<output_dir>] [--prefixes=<prefixes>] [--max-prefix-err=<max_prefix_err>] [-v | -vv | -vvv]
  freebarcodes generate     <barcode_length> <num_errors> [--output-dir=<output_dir>] [--cont=<prev_bc_fpath>] [--exclude=<exclude_bc_fpath>] [--4sets] [--cont_4sets=<prev_4sets_fpath>] [-v | -vv | -vvv]
  freebarcodes prune        <raw_barcodes_file> <num_errors> [--output-dir=<output_dir>] [-v | -vv | -vvv]
  freebarcodes concatenate  <barcode_files> [--output-dir=<output_dir>] [--max_bc=<max_bc>] [-v | -vv | -vvv]

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

Most users will be able to use one of the pre-generated lists of barcodes included in  `freebarcodes/barcodes/` for their experiments. Barcodes are stored in lists according to barcode length and number of errors corrected. For example, barcodes of length 15 which correct up to 2 errors are stored in the file `barcodes15-2.txt`. Each line of the file contains a unique barcode. The `freebarcodes/barcodes/alt_lists` folder contains alternate, non-intersecting lists of barcodes for some of those in `freebarcodes/barcodes/`.

#### Decode

Decoding barcodes is as easy as telling `freebarcodes.py` which barcode file(s) to use and which fastq file(s) should be decoded. For concatenated barcodes, list all files of barcodes in order, separated by commas and no spaces. Multiple fastq files with the same barcode list can similarly be comma-delimited. For example,
```
freebarcodes.py  decode  barcodes15-2.txt,barcodes15-2.txt  seqdata1.fastq,seqdata2.fastq
```
will decode the barcodes in `seqdata1.fastq` and `seqdata2.fastq`, where the decoding process looks for two 15 bp, 2-error correcting barcodes at the beginning of each read. 

The `--prefixes` option allows the user to look for barcodes after some known starting sequence, for example, after a primer sequence. In the event of multiple primer sequences being possible, comma-delimiting is also allowed and the observed prefix sequence is also output. Only one prefix sequence is allowed per read, however. Note that in some cases the same could be accomplished by using concatenated barcodes, but this method is more general, preferred when the possible prefixes are few (for example, 1 or 2), and possibly long and/or of different lengths. Note that the user must also specify the maximum allowed errors in each prefix with the `--max-prefix-err` option.

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

To generate lists using the above filtering criteria but excluding a certain set of barcodes, make a file with one undesired barcode per line and use the `--exclude` flag.


#### Prune

To prune a set of pre-existing barcodes to a subset of valid FREE divergence barcodes, make a file with one barcode per line and pass it `freebarcodes.py prune`. The user must also specify a number of errors to correct.

#### Concatenate

For experiments where large numbers of barcodes are needed, one can concatenate barcodes and decode them with the decode command above. Any set of barcodes can be concatenated with any other set of barcodes, but in general one will want to filter out certain barcodes that do not go well together, for example, those prone to internal hairpin structure. To concatenate any set of barcodes with any other set or sets of barcodes, filtering according to the same sequencing and synthesis filters as used for barcode generation, use the `freebarcodes.py concatenate` command with the lists of barcodes to concatenate, in order, and separated by commas and no spaces. Note that GC content is not checked, however, as it is already presumed balanced as desired. Barcodes are output as tab delimited sub-barcodes.

#### Barcode 4-sets
A friend recently requested bespoke barcodes for a microwell plate experiment he was conducting, and I share them here for their interesting properties. These barcodes come in sets of 4, here called "barcode 4-sets", where each 4-set is used in a different microplate well. The key property of each 4-set is that the barcodes must have a different base pair in each position, covering all four possibilities between the four barcodes. For example, 
```
ACTTGTTCGT             AGTAACCATG
CAAGCAGACA     and     CAAGTAATGC
GTCCTCAGTG             GCCTCTGCAA
TGGAAGCTAC             TTGCGGTGCT
```
are two valid barcode 4-sets: each column contains all four possibilities, ACGT. A valid collection of barcode 4-sets---`barcode_4sets/barcode_4sets_10-1.txt`, for example---can then all be pooled together for sequencing as a valid collection of FREE barcodes.

These barcodes are by default generated with the same sequence restrictions as above (GC content 40-60%, etc.), as well as the additional restriction that no two barcodes in the same 4-set can have more than 4 nucleotides reverse complementarity. 

Pre-generated lists can be found in the `barcode_4sets` folder, where each list has two files, one with the barcodes listed as 4-sets, and one with the barcodes listed alphabetically. The second, alphabetical file should be used for barcode decoding with the normal decoding pipeline, while the first is for experimental setup. Filenames give the barcode length and number of correctible errors. These barcodes were generated using the `freebarcodes generate` command with the `--4sets` flag.

### Examples
A couple example fastq files are included in the `examples` folder to demonstrate decoding. The first one is a "standard" example with a simple 8 bp, 2-error correcting barcode at the beginning of each sequence, given in the file `example_8-2_barcodes.fq`. Assuming the `examples` directory as the working directory, the decode command is
```
freebarcodes decode ../barcodes/barcodes8-2.txt example_8-2_barcodes.fq
```
which outputs the file `example_8-2_barcodes_decoded.txt`.

The second example, `example_8-2_barcodes_prefixes_TCTACTCTCCATACG_CACTTGGATC.fq`, shows proper usage of the prefix - detecting feature. Again, we use the length 8, 2-error correcting code, but this time we have two "primer" sequences, exactly one of which is present in each read before the barcode. Note that the two primer sequences are different lengths: 15 and 10. We also must specify a maximum number of errors to allow, for which we choose 3 and 2 respectively. The command is then
```
freebarcodes decode ../barcodes/barcodes8-2.txt example_8-2_barcodes_prefixes_TCTACTCTCCATACG_CACTTGGATC.fq --prefixes=TCTACTCTCCATACG,CACTTGGATC --max-prefix-err=3,2
```
