# make sure to compile the required file before running this script

#!/bin/bash
EXECUTABLE="./a.out" # change if needed
MAX_THREADS=4

for threads in $(seq 1 $MAX_THREADS); do
  echo -e "\nTesting with $threads threads..."
  export OMP_NUM_THREADS=$threads
  $EXECUTABLE
done
