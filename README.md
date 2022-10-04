# CNVision
R Shiny App for the visualization of SMURFSeq CNV-profiles

## Scope

This app can display the long-form output produced by the [SMURFSeq](https://github.com/smithlabcode/smurfseq_scripts.git) protocol. An implementation of this pipeline that makes use of [singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/#) features can be found [here](https://github.com/rmolbrich/SMURFsnake) as well.


## Functionality

Files are by default imported from the `data/cnv_profile` directory, located within the apps' workspace or from a user-defined space. In addition, imported data sets may be saved locally to speed up the start-up procedure for the application.

Included references comprise 'hg19' and 'hg38'. The provided data is required to enhance the plotting functionality.

## WiP

This project is a work in progress, and the implemented functionality is limited to what is sufficient for our purposes. Please raise an issue or contact me if you wish for additional functionality. 

