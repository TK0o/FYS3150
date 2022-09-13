all the .cpp files are run from terminal by first compiling and linking by 
g++ file_name.o -o file_name.exe
before then executed with
./file_name.exe
many of the .py files requires the .cpp files to be executed first, because they write the required .txt
files to the /data folder which is the directory the .py files read from
