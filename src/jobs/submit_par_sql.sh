#!/bin/sh

qsub -v dataset="1",version="par",samples="50" job_sql
qsub -v dataset="2",version="par",samples="50" job_sql
qsub -v dataset="4",version="par",samples="50" job_sql
qsub -v dataset="8",version="par",samples="50" job_sql
qsub -v dataset="16",version="par",samples="50" job_sql
qsub -v dataset="32",version="par",samples="50" job_sql
