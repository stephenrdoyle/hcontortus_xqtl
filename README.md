# hcontortus_xqtl

This repository contains the analysis workflow for the manuscript under preparation
tentatively titled "Genomic landscape of drug response reveals novel mediators of anthelmintic resistance".

Abstract:
"Like other pathogens, parasitic helminths can rapidly evolve resistance to drug treatment. The genetic basis of anthelmintic resistance is poorly understood but is key to developing methods to track and combat the spread of resistance and design new compounds that can break resistance. Using a genetic cross between drug-sensitive and drug-resistant strains of the globally important parasitic nematode Haemonchus contortus, we map resistance loci for three commonly used broad-spectrum anthelmintic classes. We identify known and reveal novel alleles for both benzimidazole and levamisole. For ivermectin, a widely used anthelmintic in animals and humans, we refine a major QTL found in resistant populations worldwide to uncover cky-1 overexpression as a potential mode of resistance. Our data define the genomic landscape of anthelmintic selection and prioritise genes and variants for the development of molecular diagnostics to combat resistance in the field.""

Authors:
**Stephen R. Doyle, Roz Laing**, David Bartley, Alison Morrison, Nancy Holroyd, Kirsty Maitland, Alistair Antonopoulos, Collette Britton, Umer Chaudhry, Ilona Flis, Sue Howell, Jennifer McIntyre, Andy Tait, Barbara Mable, Ray Kaplan, Neil Sargison, Matthew Berriman, Eileen Devaney, James A. Cotton
- Stephen and Roz are co-first authors
- Eileen and James are co-last authors


The code for the analysis is broken up into subsections as follows:  
- [Core informatics](https://github.com/stephenrdoyle/hcontortus_xqtl/blob/master/03_code/hcontortus_xqtl.workbook.mapping_variants.md)
- [Genome-wide analyses](https://github.com/stephenrdoyle/hcontortus_xqtl/blob/master/03_code/hcontortus_xqtl.workbook.genomewideplots.md)  
- [Benzimidazole analyses](https://github.com/stephenrdoyle/hcontortus_xqtl/blob/master/03_code/hcontortus_xqtl.workbook.benzimidazole.md)  
- [Levamisole analyses](https://github.com/stephenrdoyle/hcontortus_xqtl/blob/master/03_code/hcontortus_xqtl.workbook.levamisole.md)  
- [Ivermectin analyses](https://github.com/stephenrdoyle/hcontortus_xqtl/blob/master/03_code/hcontortus_xqtl.workbook.ivermectin.md)  
- [US field samples](https://github.com/stephenrdoyle/hcontortus_xqtl/blob/master/03_code/hcontortus_xqtl.workbook.US_field_genomewide.md)  

To recreate the figures, some data will need to be downloaded from my FTP repository, which can be accessed here: ftp://ftp.sanger.ac.uk/pub/project/pathogens/sd21/hcontortus_xqtl/

Note: while much effort is made to make all data and code available, some of the code is very specific to the Sanger HPC environment and will need to be modified to work elsewhere. If you are interested in using any of it for your own work, but are stuck, please do get in touch.

Any reuse of data or code is encouraged with due acknowledgement, either via citation of the published manuscript (when available) and/or GitHub repository. Comments, suggestions, and discussion are welcome.

******
## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
