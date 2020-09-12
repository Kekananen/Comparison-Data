# Rqtl
For 2019 Analysis Project with introgression regions 
The workflow of files into scipts goes as follows:
  Input: HapMap_filtered.csv, SNPs_BLASTv.4_TwoRows_DMz18-205.csv 
      -> HapMapHomoVHet.py 
        -> Output: SNPs_HAPMAP_DMz18-205v2.csv, parentsID.csv, SNPs_HapMap_DMz18-205v3.csv
        
 Input: SNPs_HapMap_DMz18-205v3.csv, Metepec_pigment_intensity.csv
      -> R:qtlConverter.py
        -> Output: rqtl.csv, Metepec_pigment_intensityv3.csv
  
