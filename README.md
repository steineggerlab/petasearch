# SRA-Search

## Installation

### Build from source

Clone this repository to your local machine. After cloning, navigate to the project folder and run the following commands.

```shell
mkdir build
cd build
cmake ..
make -j 16
```

Afterwards, the built binary `srasearch` can be found in `./build/src/`.

## Usage

### Preindexing

`srasearch` requires doing pre-indexing on the target database by calling sub-command 
`createkmertable` first. You can use the following command to generate the kmer table and ID table for a target 
database called `targetDB`.

```shell
srasearch createkmertable targetDB target_kmertable tmp
```

### Searching

```shell
srasearch petasearch queryDB targetlist resultlist tmp
```

`targetlist` should be a file containing all target databases and kmer tables. An example `targetlist` would look 
like this:

```text
targetDB1   target_kmertable1
targetDB2   target_kmertable2
```

## To-Dos

- [ ] Supports direct input one database