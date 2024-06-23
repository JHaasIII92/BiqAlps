#!/bin/bash

mpirun -np 4 ./Biq -Alps_instance /lstr/sahara/biqcrunch/TestData/MaxCut/bc/bqp50-1.sparse.bc -param Biq.par 
mpirun -np 4 ./Biq -Alps_instance /lstr/sahara/biqcrunch/TestData/MaxCut/bc/bqp50-2.sparse.bc -param Biq.par
mpirun -np 4 ./Biq -Alps_instance /lstr/sahara/biqcrunch/TestData/MaxCut/bc/bqp50-3.sparse.bc -param Biq.par
mpirun -np 4 ./Biq -Alps_instance /lstr/sahara/biqcrunch/TestData/MaxCut/bc/bqp50-4.sparse.bc -param Biq.par
mpirun -np 4 ./Biq -Alps_instance /lstr/sahara/biqcrunch/TestData/MaxCut/bc/bqp50-5.sparse.bc -param Biq.par
mpirun -np 4 ./Biq -Alps_instance /lstr/sahara/biqcrunch/TestData/MaxCut/bc/bqp50-6.sparse.bc -param Biq.par
mpirun -np 4 ./Biq -Alps_instance /lstr/sahara/biqcrunch/TestData/MaxCut/bc/bqp50-7.sparse.bc -param Biq.par
mpirun -np 4 ./Biq -Alps_instance /lstr/sahara/biqcrunch/TestData/MaxCut/bc/bqp50-8.sparse.bc -param Biq.par
mpirun -np 4 ./Biq -Alps_instance /lstr/sahara/biqcrunch/TestData/MaxCut/bc/bqp50-9.sparse.bc -param Biq.par
mpirun -np 4 ./Biq -Alps_instance /lstr/sahara/biqcrunch/TestData/MaxCut/bc/bqp50-10.sparse.bc -param Biq.par