#!/bin/sh

#qsub -v dataset="1",version="seq",samples="50" job_sql
#qsub -v dataset="2",version="seq",samples="50" job_sql
#qsub -v dataset="4",version="seq",samples="50" job_sql
#qsub -v dataset="8",version="seq",samples="50" job_sql
#qsub -v dataset="16",version="seq",samples="50" job_sql
qsub -v dataset="32",version="seq",samples="50" job_sql
