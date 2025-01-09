### simulate simple genotypes for PRS simulations ###
plink2_path=/u/home/k/kaia/STATSGEN/plink/plink2
out_path=/u/home/k/kaia/GESTALT/data/sim/genotypes/N50kM1m # directory to store simulated genotypes

### simulation parameters ###
    # N = 50,000, M = 1,000,000
    # 1 chromosome
    # phased 
    # allele frequencies are uniformly distributed in [0, 1]
    # explicit LD structure (defined by plink2)
$plink2_path --dummy 50000 1000000 acgt pheno-ct=0 phase-freq=1 \
 --make-bed \
 --out ${out_path}/sim_N50kM1m


 ### more realistic simulations ### 
    # N = 50,000, M = about 1,000,000 HAPMAP3 SNPs
    # 22 chromosomes
    # phased 
    # allele frequencies and LD structure from reference panel 1kg 
    # PRS-CSx paper and BridgePRS paper use HAPGEN2 to simulate genotypes
        # requires some work to prepare haplotypes and recombination maps because documentation is old and link to reference files are broken