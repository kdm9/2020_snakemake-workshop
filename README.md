# An introduction to workflow management


Today we're learning how to manage computational workspaces. I'll start with some broad principles around good data hygiene, which are broadly applicable to any computational analyses. We'll touch on ways to neatly organise code and data, how to version our workspaces with git, and how to take notes that remind ourselves of progress.

Then, we'll move to cover workflow management, which is a fancy way of saying "smart ways of working out which tasks are needed to complete an analysis". Specifically, we'll look at Snakemake, a python-based workflow manager with a billion fancy features that makes this task less daunting, especially for realistically-complex workflows. 

As our case-study, we'll do a very simple variant calling workflow. This consists of sequence read QC, alignment to a reference, data reformatting, and finally variant calling and filtering. We'll gloss over a lot of the details of the workflow for time reasons, but it has been very nicely outlined in [this data-carpentry tutotrial](https://gwu-omics2019.readthedocs.io/en/latest/variant_calling.html)


## Workspaces


A workspace is the directory within which some computational analyses happen. If you're used to working in R, an Rstudio project is a great example of a workspace, and of a set of tools and conventions that help us keep track of our code, data, metadata, and results.

Let's discuss key features of a well-organised workspace. Ideally, our workspaces should be:

- **Self-contained:** keep data, metadata, code, and outputs all under one project-specific directory. This also helps keep our workspaces **Sharable**.
- **Sharable:** one should be able to send (data volume notwithstanding) a workspace to a supervisor/collaborator/colleague with minimal fuss, and for them to be able to run this workflow and achieve the same results. Sharablity is a requirement for Reproducibility.
- **Reproducible:** Someone else (a reviewer, for example) should be able to convince themselves of your findings by running your code themselves, possibly with modifications. Full byte-for-byte reproducibility is less important than being able to quickly install all software, download all data, and then run the exact analysis steps that you ran.

Some recommendations:

- A directory should keep recreateable data separate from raw data. I suggest two direcrtories, `./rawdata` and `./outputs`.
- Step 0 of all analyses is to BACK UP `rawdata` NOW!!! Be a data squirrel.
- `./outputs` should ONLY contain recreateable data, such that you can do `rm -rf ./outputs` with impunity.
- `./outputs` should have subdirectories for each "step"
- Document why things are being done the way they are with comments in scripts/Snakefiles
- Use jupyter notebooks, markdown files, or labarchives as a lab notebook. Notebooks and scripts have different purposes -- notebooks are for you, scripts are for the computer.
- Use git to track code, metadata, and documentation, but not log files or data. We'll cover this later
- Use tmux to maintain persistent sessions on a remote server, one per project

## Our variant calling workflow

For the rest of today, we'll focus on a variant calling workflow. This relatively simple workflow consists of several conceptual steps, each of which are made up of one or more commands. Step 2 depends on the output of step 1, and step 3 on step 2, etc. I have already implemented this workflow for you, which you can download and extract using the following commands:

```bash
wget ... #TODO
tar xvf ...
```

These large conceptual steps are:

1. Raw sequence read data pre-processing per sample (`01_qc.sh`)
2. Alignment of reads to a reference with a short read aligner per sample (`02_align.sh`)
3. Collation of all per-sample aligned reads into one large merged BAM (`03_merge.sh`)
4. Variant calling and filtering on this merged set of samples (`04_varcall.sh`)

The details within each step are largely out of scope, but generally cobble together a series of smaller computational tasks that together perform one of these larger conceptual steps.


# Snakemake

For the remainder of today, we'll work on translating our variant-calling workflow above to a snakemake workflow. We'll follow the excellent online tutorial on snakemake, which you can find at <https://snakemake.readthedocs.io/en/stable/tutorial/basics.html>.
