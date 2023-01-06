#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <sys/mman.h>

int main (int argc, char *argv[]){
    FILE* handle = fopen(argv[1], "rb");
    int fd=fileno(handle);
    struct stat fileStat;
    fstat(fd, &fileStat);
    size_t fileSize = fileStat.st_size;

    long* pos = (long*)mmap(NULL, fileStat.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    // madvice MADV_SEQUENTIAL
    if (posix_madvise (pos, fileStat.st_size, POSIX_MADV_SEQUENTIAL|POSIX_MADV_WILLNEED) != 0){
        std::cout << "posix_madvise returned an error " << "\n";
    }
    
    size_t count=0;
    long sum = 0;
    long lastRead = 0;
    for(size_t i = 0; i < fileSize/sizeof(long); ++i){
        count++;
        lastRead = *(pos+i);
        sum += lastRead;
    }
    std::cout<<"filesize/8: "<< fileSize/8<<std::endl;
    std::cout<<"lastRead: "<<lastRead<<std::endl;
    std::cout<<"sum:" <<sum<<std::endl<<"count: "<<count<<std::endl;
    
    munmap(pos,fileSize);
    fclose(handle);

    count = 0;
    lastRead = 0;
    sum = 0;
    std::ifstream myFile ("../test_k-merTable", std::ios::in | std::ios::binary);
    for(size_t i = 0; i < fileSize/sizeof(long); ++i){
       count ++;
       myFile.read((char *) &lastRead,sizeof(long));
       sum+= lastRead; 
    }
    myFile.close();
    std::cout<<"sum:" <<sum<<std::endl;
    std::cout<<"lastRead:" <<lastRead<<std::endl<<"count: "<<count<<std::endl;
}
