Workflow for ntuple generation

`run_module.py` : test file for latest configuration

`misc/condorJobMaker.py` : can produce condor jobs for making ntuples
`misc/makeJobs.sh`       : script that uses `misc/condorJobMaker.py` to make the jobs .. modify it to your need 
`misc/genNtuplize.py`    : the config file used by `misc/makeJobs.sh`
`fileList/c3_1_c4_1_HHHto4b2gamma.fls` : file list used to make the jobs


To produce the condor jobs
```bash
./misc/makeJobs.sh
ls Condor/Job*/*.sub
#Identify ur submit file
condor_submit <ur file>
```
