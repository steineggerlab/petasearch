# Petasearch

Petasearch enables searching through the largest sets of proteins.

## Installation

`Petasearch` depends on [`block-aligner`](https://github.com/Daniel-Liu-c0deb0t/block-aligner) for fast computation
of Smith-Waterman alignments in the `blockalign` module. Thus, [the Rust Programming 
Langugage](https://www.rust-lang.org/) needs be installed on the user's machine.

### Build from source

Clone this repository to your local machine. After cloning, navigate to the project folder and run the following commands.

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 16
```

Afterwards, the built binary `srasearch` can be found in `./build/src/`. You can add the path to the binary to your `$PATH`
in order to user it easily.

## Usage

### Create Petasearch databases

To achieve space efficiency, `srasearch` will store the target databases in a specific highly-compressed format. You can use `convert2sradb` to convert a FASTA/FASTQ file or a MMseqs2 database into a srasearch database.

Example usage:

```shell
srasearch convert2sradb target.fasta targetDB
srasearch convert2sradb mmseqsDB targetDB
```

### Preindexing

`srasearch` requires pre-indexing the target database by calling the module 
`createkmertable` first. You can use the following command to generate the kmer table and ID table for a target
database called `targetDB`.

```shell
srasearch createkmertable targetDB target_kmertable
```

### Combined workflow

`Petasearch` provides a combined workflow that will produce only one output file `alignments.m8` that contain all the 
search results from searching `queryDB` against all the target databases listed in `targetlist`.

```shell
srasearch petasearch queryDB targetlist resultlist alignments.m8 tmp
```
### Easy workflow

`Petasearch` also provides an easy workflow that will accept `fasta` file as the input query dataset. The user also do
not need to provide `resultlist` as an input.

```shell
sraserach easy-petasearch query.fasta targetlist alignments.m8 tmp
```

### Use index table to prefilter query-target pairs

To search a MMseqs2 database against a list of Petasearch databases, simply run:

```shell
srasearch comparekmertables queryDB targetlist resultlist
```

`targetlist` should be a file containing all target databases and kmer tables. An example `targetlist` would look
like this:

```text
target_kmertable1   targetDB1
target_kmertable2   targetDB2
```

`resultlist` should be a file containing all file names for output files to store the prefiltering result. An example
`resultlist` would look like this:

```text
compkmer_res_1
compkmer_res_2
```

### Compute Smith-Waterman alignment selectively

```shell
srasearch blockalign queryDB targetDB1 compkmer_res_1 compali_res_1
```

### Print out alignment results

```shell
srasearch convertsraalis queryDB targetDB1 compali_res_1 alignments_1.m8
```

`resultlist` should be a file with the *same* amount of entries of `targetlist`. Those are the file names

