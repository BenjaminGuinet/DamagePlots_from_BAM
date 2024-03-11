# BAM TO DAMAGE PLOTS 

## Description

This script converts a bam file into a PDF file with multiple plots including reads as well as damage informations.

You will need Amber, Picard, samtools, PDFtools as well as multiple python packages (see header of the .py script) to be installed. 

## Usage Example

```
python3 Make_multiple_plotV2.py -b Yersinia_pestis_control.bam 
```
______________

## Main outputfiles :

- ***NoDup_sorted_filterBAM_quality_filtred_Yersinia_pestis_control_dna_damage.csv*** :  CSV file output with many important informations.
- ***NoDup_sorted_filterBAM_quality_filtred_Yersinia_pestis_control_dna_damage.pdf*** :  PDF file output with the plots.

Other files correspond to amber and pydamage outputs. 
______________

## Output PDF example : 

[View the PDF](https://github.com/BenjaminGuinet/DamagePlots_from_BAM/blob/main/NoDup_sorted_filterBAM_quality_filtred_Yersinia_pestis_control.bam_dna_damage_plots.pdf)



