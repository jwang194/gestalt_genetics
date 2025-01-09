scripts to simulate phenotypes and gwas 

loop_configs.sh 
    loops through every possible parameter combinations
        for now all 10 traits 10% heritability
    polygenicity: [0.01, 0.05, 0.1]
	proportion of shared causal variants [0.05, 0.2, 0.4]
	genetic covariance between shared components of traits [0.01, 0.04, 0.08]
	environmental covariance between traits [-0.1 0.1] 

Simulation design: 
    polygenicity controls the number of total causal variants
    proportion of shared causal variants extracts a fraction of the causal variants to be shared between traits
        the other causal variants are specific to a trait, partitioned equally across traits 
    the additive genetic variance of the shared and specific components are set to be the same value (heritability)
        the betas are scaled so that we get the heritability we want 
        so the total genetic covariance between traits is half of what we set it to
    in this simulation, the per-SNP heritability is determined by the number of causal variants that are shared or specific

    here, the maximum total genetic covariance we can simulate is 0.5 in this setup
        what does it mean for traits to have genetic correlation 0.8? 
            genetic covariance is 80% of the sqrt(genetic variance trait 1)*sqrt(genetic variance trait 2)
            means a large majority of causal variants are shared (with high correlation between effect sizes) 
            if there are few shared causal variants, but high genetic correlation, the specific effects sizes must be small?

        for example:
            2 traits
            trait heritability = additive genetic variance = 0.1  (assume total phenotypic variance = 1)
            genetic correlation = 0.8, genetic covariance = 0.08
            1,000 total causal variants, 200 shared, 400 specific to each trait 

            in an extreme case, if we assume the effect sizes for the shared variants are identical
                genetic correlation of shared component = 1
                genetic covariance of shared component = heritability = 0.08
            then the specific component of each phenotype should have the following genetic variance: 
                0.8 = 0.08/(0.08 + x)
                x = (0.08 - 0.08*0.1)/0.08 = 0.02
                the trait-specific variant effect sizes are drawn from mean 0 variance 0.02/400 (# trait-specific snps)
                    which is small compared to shared variant effect sizes of variance 0.08/200 (# shared snps)


    need to understand these assumptions to improve this simulation framework. 

    move with what we have for now, come back

Better design? 
    keep our current design 
        the input genetic covariance matrix determines the genetic covariance of the SHARED components 
    add another parameter to dictate the total genetic covariance between traits 
        this determines the additive genetic effect of the trait-specific component
            and therefore the per-snp heritability of specific vs shared effects
        the total genetic covariance input cannot be greater than the genetic covariance of the SHARED component
            if we want total genetic covariance = genetic covariance of the SHARED component, this means additive genetic variance of trait-specific component is 0

    OR 

    input shared and specific fraction, make the assumption that all causal variants have the same per-SNP heritability 
    the additive genetic variance of the shared and specific components now depend on the number of SNPs in each category 
    input the desired genetic correlation 
        set genetic covariance of shared effects to r_g * sqrt(genetic variance trait 1)*sqrt(genetic variance trait 2)
        the genetic covariance cannot be greater than the additive genetic variance of the shared component 
            so genetic correlation is at most (# shared variants / # total variants)

    
    
    
