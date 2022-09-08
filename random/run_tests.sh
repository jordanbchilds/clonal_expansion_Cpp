#!/usr/bin/env bash


echo -e "\n---- DUMMY --------------------------------------------\n"

make clean && make
echo
./random_stream

echo -e "\n---- URAND32 ------------------------------------------\n"

make clean && make EXTRA_FLAGS=-DUSE_URAND32=1

echo -e "\nWith non-random PRNG state (0,0,0,0)\n"

./random_stream 0 0 0

echo -e "\nWith random PRNG state\n"

./random_stream 0

echo -e "\n---- URAND_F32 ----------------------------------------\n"

make clean && make EXTRA_FLAGS=-DUSE_URAND_F32=1

echo -e "\nWith non-random PRNG state (0,0,0,0)\n"

./random_stream 0 0 0

echo -e "\nWith random PRNG state\n"

./random_stream 0

echo -e "\ndone\n"

