# SRA-Search

## How to build the project

After cloning, navigate to the project folder and run the following commands to use `CMake` to automate the building process.

```shell script
mkdir build
cd build
cmake ..
make -j 16
```

This project utilizes quite a few features provided by `mmseqs2`. So also remember to clone [MMseqs2 on Github](https://github.com/soedinglab/MMseqs2.git). And then navigate into the target folder and run the commands listed above.

Afterwards, we should build the database needed for testing.

```shell script
mkdir /path/to/your/database
cd /path/to/mmseqs2/build/
src/mmseqs databases UniProtKB/Swiss-Prot path/to/your/database/swissprot tmp
```
And then we can run `createKmerTable` to test that 

```shell script
cd /path/to/SRA-search/build
src/srasearch createkmertable /path/to/your/database/swissprot path/to/your/database/swissprot_kmertable -k 8
```

## Goal

While running 
```shell script
src/srasearch createkmertable <database> <database_kmertable> -k 11
```
We can clearly feel that:

1. The kmerTable generated is taking a hugh amount of disk space (more than 6GB)
2. The sorting algorithm is very slow

We would like to optimize the space complexity by making the kmerTable created taking less disk space, also we would like to speed up the sorting period. 

While using CLion to build the project, set the program arguments of `sraserach|Debug` as:

```shell script
createkmertable /home/matchy233/summer-research/database/swissprot /home/matchy233/summer-research/database/swissprot_kmertable -k 8
```

If after the modifications the program still runs nicely under the debug condition, we should change `k` to `11` and the see if there's any improvement in the disk space usage!