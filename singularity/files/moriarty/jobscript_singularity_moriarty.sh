#!/bin/sh
# properties = {properties}
singularity exec  -B /usr/bin/:/moriarty/bin/ -B /gluster-storage-volume,/data,/etc/passwd,/sw,/run/munge,/usr/lib64,/projects caps.simg  bash -c "{exec_job}"
