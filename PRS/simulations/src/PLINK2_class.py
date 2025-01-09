import subprocess
import pandas as pd

### python class for plink2 object ### 
    # stores paths/inputs for flags as parameters
    # functions represent various plink2 operations, automatically using parameters as inputs to flags 
        # do not need to manually specify inputs paths for each flag 
    # can update the parameters 
        # useful when processing genotype files in multiple steps 

class Plink2:
    def __init__(self, plink2_path, bfile, threads="1"):
        self.plink2_path = plink2_path
        self.bfile = bfile
        self.threads = str(threads)

    def set_param(self, name, value):
        # can be used to add additional parameters to plink2 object (ex. phenotype file)
        setattr(self, name, value)

    def simulate_genotypes(self, N, M, output):
        # N = number of individuals, M = number of SNPs
        # simulates phased genotypes with allele frequencies drawn from uniform distribution [0, 1]
        # plink2 documentation describes that explicit LD is simulated 
        cmd = [self.plink2_path, 
        "--dummy", str(N), str(M), "acgt", "pheno-ct=0", "phase-freq=1", 
        "--threads", self.threads,
        "--out", output]
        subprocess.run(cmd, text=True)
        ### creates and returns new plink2 object with updated bfile path ###
        return Plink2(self.plink2_path, output, self.threads)
        
    def keep_ids(self, ids, output):
        # ids: file path, first two columns FID/IID
        cmd = [self.plink2_path, 
        "--bfile", self.bfile, 
        "--threads", self.threads,
        "--make-bed", 
        "--keep", ids,
        "--out", output]
        subprocess.run(cmd, text=True)
        ### creates and returns new plink2 object with updated bfile path ###
        return Plink2(self.plink2_path, output, self.threads)

    def remove_ids(self, ids, output):
        # ids: file path, first two columns FID/IID
        cmd = [self.plink2_path, 
        "--bfile", self.bfile, 
        "--threads", self.threads,
        "--make-bed", 
        "--remove", ids,
        "--out", output]
        subprocess.run(cmd, text=True)
        ### creates and returns new plink2 object with updated bfile path ###
        return Plink2(self.plink2_path, output, self.threads)

    def filter_phenotype(self, P_df):
        # P_df: dataframe of phenotypes first two columns FID/IID

        # This function filters the phenotype to only keep individuals in self.bfile + '.fam' file
        # Normalizes the phenotypes and saves the file path to self.pheno flag 
        fam = pd.read_csv(self.bfile + '.fam', sep = '\t', header = None)
        P_df_filtered = P_df[P_df.iloc[:,1].isin(fam.iloc[:,1])]
        P_df_filtered = P_df_filtered.set_index(P_df_filtered.columns[1], drop = False)
        P_df_filtered = P_df_filtered.reindex(fam.iloc[:, 1]).reset_index(drop = True)
        scaled_P = P_df_filtered.copy()
        scaled_P.iloc[:, 2:] = (P_df_filtered.iloc[:, 2:] - P_df_filtered.iloc[:, 2:].mean()) / P_df_filtered.iloc[:, 2:].std()
        scaled_P.to_csv(self.bfile + '.pheno', sep = '\t', index = False)
        self.set_param('pheno', self.bfile + '.pheno')

    def compute_freq(self):
        # computes allele frequencies and saves to self.bfile + '.acount'
        cmd = [self.plink2_path, 
        "--bfile", self.bfile, 
        "--threads", self.threads,
        "--freq", "counts",
        "--out", self.bfile]
        subprocess.run(cmd, text=True)

    def gwas(self, freq): # currently only expects continuous phenotypes
        # gwas on pheno file 
            # use self.filter_phenotype above to ensure filtered phenotype file is saved to self.pheno
        # returns an dictionary of phenotype name <-> results file paths, an entry for each phenotype
        cmd = [self.plink2_path, 
        "--bfile", self.bfile, 
        "--threads", self.threads,
        "--pheno", self.pheno,
        "--read-freq", freq,
        "--glm", "allow-no-covars",
        "--no-fam-pheno",
        "--out", self.bfile]
        subprocess.run(cmd, text=True)

        # collect results files 
        results_map = {}
        P_df = pd.read_csv(self.pheno, sep = '\t')
        for name in P_df.columns[2:]:
            results_map[name] = '%s.%s.glm.linear'%(self.bfile, name)

        return results_map


    def clump(self, gwas_file, snp_field, p_field, out):
        # function to use this plink2 object as the LD reference to clump/prune association study results
        # uses default parameters used in PRSice-2 and PRS-CSx paper
        cmd = [self.plink2_path, 
        "--bfile", self.bfile, 
        "--threads", self.threads,
        "--clump", gwas_file,
        "--clump-p1", "1",
        "--clump-r2", "0.1",
        "--clump-kb", "250",
        "--clump-id-field", snp_field,
        "--clump-p-field", p_field,
        "--out", out]
        subprocess.run(cmd, text=True)
        cmd = "awk 'NR!=1{print $3}' %s.clumps > %s.valid.snps"%(out, out)
        subprocess.run(cmd, text=True, shell=True)
        return "%s.valid.snps"%(out)

    def create_score_range(self, thresholds):
        # thresholds = array of p-value thresholds
        # creates a valid file for --q-score-range and saves to score_range parameter
        range_list = []
        for threshold in thresholds:
            range_list.append('%s 0 %s'%(str(threshold), str(threshold)))
        pd.DataFrame(range_list).to_csv( self.bfile + '.score_range', index = False, header = None )
        self.set_param('score_range', self.bfile + '.score_range' )

        #assoc_df = pd.read_csv(assoc_file, sep = '\t')
        #assoc_df[[snp_field, p_field]].to_csv(self.bfile + '.score_range_data', index = False, sep = '\t')
        #self.set_param('score_range_data', self.bfile + '.score_range_data' )
        


    def prs(self, gwas_file, snp_field, p_field, a_field, B_field, clumps_file, p_thresholds, out):
        ### INPUTS ###
            # gwas_file = path to gwas file 
            # snp_field = name of snp id column in gwas file
            # p_field = name of p-value column in gwas file 
            # a_field = name of effect allele column in gwas file 
            # B_field = name of effect size column in gwas file 
            # clumps_file = path to clumped snp ids 
            # p_thresholds = array of p-value thresholds 
            # out = path+prefix to output file

        ### OUTPUT ###
            # prs_map = dictionary mapping p-value threshold <-> dataframe of prs scores

        # create q-score-range file and extract snp id and p-value for q-score-range data file
        self.create_score_range(p_thresholds)

        # extract column numbers for each field 
        gwas_header = pd.read_csv(gwas_file, sep = '\t', nrows = 0)
        snp_index = str(gwas_header.columns.get_loc(snp_field) + 1)
        allele_index = str(gwas_header.columns.get_loc(a_field) + 1)
        beta_index = str(gwas_header.columns.get_loc(B_field) + 1)
        p_index = str(gwas_header.columns.get_loc(p_field) + 1)

        # prs scoring 
        cmd = [self.plink2_path, 
        "--bfile", self.bfile, 
        "--threads", self.threads,
        "--score", gwas_file, snp_index, allele_index, beta_index, 'header',
        "--q-score-range", self.score_range, gwas_file, snp_index, p_index, 'header',
        "--extract", clumps_file,
        "--out", out]
        subprocess.run(cmd, text=True)

        prs_map = {}
        for value in p_thresholds:
            scores_df = pd.read_csv('%s.%s.sscore'%(out, value), sep = '\t')
            #prs_map[value] = scores_df.iloc[:, [0,1,5]]
            prs_map[value] = scores_df.iloc[:, [0,1,-1]]

        return prs_map

