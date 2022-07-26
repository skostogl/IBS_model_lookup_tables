### Instructions from G. Sterbini from https://github.com/xsuite/example_DA_study/tree/release/v0.1.0 and IBS scripts from F. Antoniou
### Compatible with tree_maker v0.1.0, submission to htcondor, lsf and slurm supported.

### Installation instructions

```bash 
# install miniconda, then:
pip install xsuite
pip install cpymad
git clone https://github.com/xsuite/tree_maker.git
pip install -e tree_maker/
git clone https://github.com/lhcopt/lhcmask.git
pip install -e lhcmask/
git clone https://github.com/lhcopt/lhctoolkit.git
git clone https://github.com/lhcopt/lhcerrors.git
```

### Connecting to the cnaf.infn.it
Go to the the hpc-201-11-01-a machine.
I do it by 
```bash
ssh bologna
```

but you need to configure your `~/.ssh/config` by adding
```bash
Host bastion
 HostName bastion.cnaf.infn.it
 User sterbini
 ForwardX11 yes

Host bologna
 ProxyCommand ssh -q bastion nc hpc-201-11-01-a 22
 ForwardX11 yes
```

As you can see, `bologna` (hpc-201-11-01-a) host is passing via the `bastion` (bastion.cnaf.infn.it) connection.


# Documentation of master study

000_make_sequence_tree.py: creates tree for all sequences from Run 3 repo. Run 002_chronjob with tree_maker_sequences.json

001_make_ibs_tree.py: creates tree for IBS scan. At the moment, 50 points of the scan per job. Run 002_chronjob with tree_maker.json

002_chronjob.py: Submits jobs in lsf, slurm or htcondor

003_save_lookup_table.py: reads output of all jobs and saves lookup table

004_interpolate.py: performs interpolation and even needed saves interpolation in pickles (I have not found a better way at the moment)

004a_sanity_check.py: compares emittance evolution derived from lookup tables with MADX 
results (in "sanity_checks" folder)

