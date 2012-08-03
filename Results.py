__author__ = 'chris'

class Results:
    def __init__(self, num_samples, num_cols, dist):
        self.mb_orig_diff = None
        self.mb_recon_diff = None
        self.paralin_cluster_diff = None
        self.paralin_nj_diff = None
        self.bc_cluster_diff = None
        self.bc_nj_diff = None
        self.r_uni_cluster_diff = None
        self.u_uni_cluster_diff = None
        self.w_uni_cluster_diff = None
        self.bc_pcoa_diff = None
        self.paralin_pcoa_diff = None
        self.num_samples = num_samples
        self.num_cols = num_cols
        self.dist = dist

    def print_results(self, out_file):
        out_file.write("%d,%d,%s,"
                       "%d,%d,%.2f,"
                       "%d,%d,%.2f,"
                       "%d,%d,%.2f,"
                       "%d,%d,%.2f,"
                       "%d,%d,%.2f,"
                       "%d,%d,%.2f,"
                       "%d,%d,%.2f,"
                       "%d,%d,%.2f,"
                       "%d,%d,%.2f,"
                       "%d,%d,%.2f\n" %
                       (self.num_samples,
                        self.num_cols,
                        self.dist,
                        self.mb_orig_diff[0],
                        self.mb_orig_diff[1],
                        self.mb_orig_diff[2],
                        self.mb_recon_diff[0],
                        self.mb_recon_diff[1],
                        self.mb_recon_diff[2],
                        self.u_uni_cluster_diff[0],
                        self.u_uni_cluster_diff[1],
                        self.u_uni_cluster_diff[2],
                        self.w_uni_cluster_diff[0],
                        self.w_uni_cluster_diff[1],
                        self.w_uni_cluster_diff[2],
                        self.paralin_cluster_diff[0],
                        self.paralin_cluster_diff[1],
                        self.paralin_cluster_diff[2],
                        self.bc_cluster_diff[0],
                        self.bc_cluster_diff[1],
                        self.bc_cluster_diff[2],
                        self.paralin_nj_diff[0],
                        self.paralin_nj_diff[1],
                        self.paralin_nj_diff[2],
                        self.bc_nj_diff[0],
                        self.bc_nj_diff[1],
                        self.bc_nj_diff[2],
                        self.bc_pcoa_diff[0],
                        self.bc_pcoa_diff[1],
                        self.bc_pcoa_diff[2],
                        self.paralin_pcoa_diff[0],
                        self.paralin_pcoa_diff[1],
                        self.paralin_pcoa_diff[2]
                        ))
        out_file.flush()







