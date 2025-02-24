mkdir -p obj
mkdir -p bin

gcc -c util.c -o ./obj/util.o
gcc -c mc.c -o ./obj/mc.o
gcc -c argp.c -o ./obj/argp.o
gcc -c globals.c -o ./obj/globals.o
gcc -c printer.c -o ./obj/printer.o

gcc -c main.c -o ./obj/main.o

gcc ./obj/*.o -o ./bin/main

# ./bin/main -s 1000 -t 0.5 -T 0.000001 -x 2 -X 0.0 -s 1024
