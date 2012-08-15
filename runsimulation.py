__author__ = 'chris'

import os
from Results import *
import app
import datetime
import stopwatch
import sys
import traceback

sys.setrecursionlimit(1000000)
timer = stopwatch.Timer()

hostname = os.uname()[1]

num_runs = 10

bits = 3

col_start = 100
col_end = 101
col_inc = 1
col_range = range(col_start, col_end, col_inc)

sample_start = 10
sample_end = 11
sample_inc = 1
sample_range = range(sample_start, sample_end, sample_inc)

dists = ['gamma', 'normal']
gamma_shape = 1
gamma_scale = 1000
mean = 10000
sd = 5000

n_gen = 100000
mpi = "/opt/local/bin/mpirun"
mb = "/Users/chris/src/mrbayes_3.2.1/src/mb"
procs = 4
project_dir = "/Users/chris/projects/bayessim"

if 'godel' in hostname:
    mpi = '/test/riveralab/cfriedline/bin/mpirun'
    mb = '/test/riveralab/cfriedline/src/mrbayes_3.2.1/src/mb'
    procs = 8
    project_dir = '/test/riveralab/cfriedline/projects/bsim2'
elif 'phylogeny' in hostname:
    mpi = '/usr/local/bin/mpirun'
    mb = '/home/cfriedline/src/mrbayes_3.2.1/src/mb'
    procs = 8
    project_dir = '/home/cfriedline/projects/bsim'

if not 'local' in hostname:
    col_start = 600
    col_end = 12001
    col_inc = 600
    col_range = range(col_start, col_end, col_inc)
    num_runs = 10

run_dir_name = datetime.datetime.now().strftime("%m%d%y_%H%M%S")
result_dir = app.create_dir(os.path.join(project_dir, "results"))
run_dir = app.create_dir(os.path.join(result_dir, run_dir_name))

out_dir = app.create_dir(os.path.join(run_dir, "out"))
log_dir = app.create_dir(os.path.join(run_dir, "log"))
log_file = open(os.path.join(log_dir, "log.txt"), "w")
app.log_file = log_file

out_file = open(os.path.join(out_dir, "out.txt"), "w")
out_file.write("num_samples,num_cols,dist,"
               "mb_orig_topo,mb_orig_symm,mb_orig_path,"
               "mb_recon_topo,mb_recon_symm,mb_recon_path,"
               "u_uni_cluster_topo,u_uni_cluster_symm,u_uni_cluster_path,"
               "w_uni_cluster_topo,w_uni_cluster_symm,w_uni_cluster_path,"
               "paralin_cluster_topo,paralin_cluster_symm,paralin_cluster_path,"
               "bc_cluster_topo,bc_cluster_symm,bc_cluster_path,"
               "u_uni_nj_topo,u_uni_nj_symm,u_uni_nj_path,"
               "w_uni_nj_topo,w_uni_nj_symm,w_uni_nj_path,"
               "paralin_nj_topo,paralin_nj_symm,paralin_nj_path,"
               "bc_nj_topo,bc_nj_symm,bc_nj_path,"
               "u_uni_pcoa_topo,u_uni_pcoa_symm,u_uni_pcoa_path,"
               "w_uni_pcoa_topo,w_uni_pcoa_symm,w_uni_pcoa_path,"
               "bc_pcoa_topo,bc_pcoa_symm,bc_pcoa_path,"
               "paralin_pcoa_topo,paralin_pcoa_symm,paralin_pcoa_path"
               "\n")
print "Running %s runs each for %s samples and %s cols" % (num_runs, sample_range, col_range)

app.create_R(out_dir)
all_results = []
completed = 0
smallest_max_otu_value = app.compute_smallest_max()
for i in range(len(col_range)):
    num_cols = col_range[i]
    taxa_tree = app.create_tree(app.find_usable_length(num_cols, bits) / bits, type = 'T')
    assert app.is_binary_tree(taxa_tree) == True
    for j in range(len(sample_range)):
        num_samples = sample_range[j]
        sample_trees = [None] * num_runs
        for dist in dists:
            results = Results(num_samples, num_cols, dist)
            all_results.append(results)
            ranges = None
            if dist == 'normal':
                ranges = app.get_range_from_normal(num_cols, bits, mean, sd, smallest_max_otu_value)
            else:
                ranges = app.get_range_from_gamma(num_cols, bits, gamma_shape, gamma_scale, smallest_max_otu_value)
            for k in range(num_runs):
                status = "%d/%d: running %d samples and %d cols (%s) %s" % (
                    completed + 1, num_runs * len(dists) * len(sample_range) * len(col_range), num_samples, num_cols,
                    dist, str(timer))

                print status
                log_file.write("%s\n" % status)
                log_file.flush()

                sample_tree = None

                #create the sample tree, if not already there and save for next distribution
                if sample_trees[k] is None:
                    sample_tree = app.create_tree(num_samples, type = 'S')
                    assert app.is_binary_tree(sample_tree) == True
                    sample_trees[k] = sample_tree
                else:
                    sample_tree = sample_trees[k]

                #setup all the matrices
                matrix = sample_names = None
                try:
                    sample_tree2, matrix, sample_names = app.create_discrete_matrix(num_cols,
                                                                                    num_samples, sample_tree, bits)
                    if sample_tree2 != sample_tree:
                        assert app.is_binary_tree(sample_tree2) == True
                        sample_trees[k] = sample_tree2
                        sample_tree = sample_tree2
                except Exception, err:
                    sys.stderr.write('ERROR! %s\n' % str(err))
                    traceback.print_exc()
                    exit(1)

                gap = app.get_range_standardized_matrix_from_discrete(matrix, bits, num_cols)
                abund = app.get_abundance_matrix(gap, ranges, dist)
                gap2 = app.restandardize_matrix(abund, ranges)
                matrix2 = app.get_discrete_matrix_from_standardized(gap2, bits, sample_names)
                matrix_cor = app.correlate_matrices(matrix, matrix2)


                #do unifrac in python (PyCogent)
                (u_uni_matrix, u_uni_rownames), (w_uni_matrix, w_uni_rownames) = app.calculate_unifrac(abund, sample_names,
                                                                                                   taxa_tree)

                # do mrbayes
                mb_tree = app.run_mrbayes(k, matrix, sample_names, num_cols, n_gen, mpi, mb, procs, dist, run_dir,
                                          num_samples, "orig")
                assert app.is_binary_tree(mb_tree) == True
                mb_tree2 = app.run_mrbayes(k, matrix2, sample_names, num_cols, n_gen, mpi, mb, procs, dist, run_dir,
                                           num_samples, "recon")
                assert app.is_binary_tree(mb_tree2) == True

                # do average linkage
                paralin_cluster_tree = app.get_paralinear_cluster()
                assert app.is_binary_tree(paralin_cluster_tree) == True
                bc_cluster_tree = app.get_bc_cluster(abund, sample_names)
                assert app.is_binary_tree(bc_cluster_tree) == True
                #unweighted
                u_unifrac_cluster_tree = app.get_py_unifrac_cluster(u_uni_matrix, u_uni_rownames)
                assert app.is_binary_tree(u_unifrac_cluster_tree) == True
                #weighted
                w_unifrac_cluster_tree = app.get_py_unifrac_cluster(w_uni_matrix, w_uni_rownames)
                assert app.is_binary_tree(w_unifrac_cluster_tree) == True

                # do neighbor joining
                bc_nj_tree = app.get_bc_nj()
                assert app.is_binary_tree(bc_nj_tree) == True
                paralin_nj_tree = app.get_paralinear_nj()
                assert app.is_binary_tree(paralin_nj_tree) == True
                u_unifrac_nj_tree = app.get_unifrac_nj(u_uni_matrix, u_uni_rownames)
                assert app.is_binary_tree(u_unifrac_nj_tree) == True
                w_unifrac_nj_tree = app.get_unifrac_nj(w_uni_matrix, w_uni_rownames)
                assert app.is_binary_tree(w_unifrac_nj_tree) == True

                # do pcoa
                bc_pcoa_tree = app.get_bc_pcoa_tree()
                assert app.is_binary_tree(bc_pcoa_tree) == True
                paralin_pcoa_tree = app.get_paralin_pcoa_tree()
                assert app.is_binary_tree(paralin_pcoa_tree) == True
                u_unifrac_pcoa_tree = app.get_unifrac_pcoa_tree(u_uni_matrix, u_uni_rownames)
                if u_unifrac_pcoa_tree is not None:
                    assert app.is_binary_tree(u_unifrac_pcoa_tree)

                w_unifrac_pcoa_tree = app.get_unifrac_pcoa_tree(w_uni_matrix, w_uni_rownames)
                if w_unifrac_pcoa_tree is not None:
                    assert app.is_binary_tree(w_unifrac_pcoa_tree)

                # calc differences
                results.mb_orig_diff = app.calculate_differences_r(sample_tree, mb_tree)
                results.mb_recon_diff = app.calculate_differences_r(sample_tree, mb_tree2)

                results.u_uni_cluster_diff = app.calculate_differences_r(sample_tree, u_unifrac_cluster_tree)
                results.u_uni_pcoa_diff = app.calculate_differences_r(sample_tree, u_unifrac_pcoa_tree)
                results.u_uni_nj_diff = app.calculate_differences_r(sample_tree, u_unifrac_nj_tree)

                results.w_uni_cluster_diff = app.calculate_differences_r(sample_tree, w_unifrac_cluster_tree)
                results.w_uni_pcoa_diff = app.calculate_differences_r(sample_tree, w_unifrac_pcoa_tree)
                results.w_uni_nj_diff = app.calculate_differences_r(sample_tree, w_unifrac_nj_tree)

                results.paralin_nj_diff = app.calculate_differences_r(sample_tree, paralin_nj_tree)
                results.paralin_cluster_diff = app.calculate_differences_r(sample_tree, paralin_cluster_tree)
                results.paralin_pcoa_diff = app.calculate_differences_r(sample_tree, paralin_pcoa_tree)

                results.bc_cluster_diff = app.calculate_differences_r(sample_tree, bc_cluster_tree)
                results.bc_nj_diff = app.calculate_differences_r(sample_tree, bc_nj_tree)
                results.bc_pcoa_diff = app.calculate_differences_r(sample_tree, bc_pcoa_tree)

                #                results.paralin_cluster_diff = (0, 0, 0)
                #                results.paralin_nj_diff = (0, 0, 0)
                #                results.paralin_pcoa_diff = (0, 0, 0)

                # print all the stuff
                app.print_matrices(abund, ranges, gap, gap2, matrix, matrix2, log_dir, k, dist, num_samples, num_cols)
                results.print_results(out_file)
                out1 = open(os.path.join(log_dir, "origtree_%d_%d_%d.tre" % (num_samples, num_cols, k)), "w")
                out2 = open(os.path.join(log_dir, "mbtree_orig_%d_%d_%s_%d.tre" % (num_samples, num_cols, dist, k)),
                            "w")
                out3 = open(os.path.join(log_dir, "mbtree_recon_%d_%d_%s_%d.tre" % (num_samples, num_cols, dist, k)),
                            "w")

                #out3 = open(os.path.join(log_dir, "unitree_%d_%d_%s_%d.tre" % (num_samples, num_cols, dist, k)), "w")
                writers = [out1, out2, out3]
                out1.write(sample_tree.as_newick_string() + ";\n")
                out2.write(mb_tree.as_newick_string() + ";\n")
                out3.write(mb_tree2.as_newick_string() + ";\n")

                #out3.write(r_unifrac_cluster_tree.as_newick_string() + ";\n")
                app.print_trees_to_pdf(taxa_tree, sample_tree, mb_tree, mb_tree2,
                                       u_unifrac_cluster_tree, w_unifrac_cluster_tree,
                                       results.mb_orig_diff, results.mb_recon_diff, results.u_uni_cluster_diff,
                                       results.w_uni_cluster_diff,
                                       dist, k)
                [i.close() for i in writers]
                completed += 1
            app.log("\t%d %s iters done in %s" % (num_runs, dist, str(timer)), log_file)
        app.log("\t%d samples done in %s" % (num_samples, str(timer)), log_file)
    app.log("\t%d cols in %s" % (num_cols, str(timer)), log_file)
app.close_R()
timer.stop()
log_file.close()
out_file.close()
print "Done in %s " % str(timer)