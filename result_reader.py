__author__ = 'chris'

import os

def main():
    #   read_results("/Users/chris/PycharmProjects/BayesianSimulation/bayesiansimulation/072412_130800/out/out.txt")
#    read_results("/Users/chris/Desktop/out.txt")

#    read_results("/Users/chris/projects/081512_105614/out/out.txt")
    read_results("/Users/chris/projects/bayessim/results/082212_104332/out/out.txt")
def read_results(file_name):
    file = os.path.abspath(file_name)
    fh = open(file)
    header = fh.readline().rstrip(os.linesep)
    results = []
    for line in fh:
        results.append(line.rstrip(os.linesep).split(","))

    process_results(header, results, os.path.dirname(file))


def process_results(header, results, out_dir):
    types = ["symm", "path", "topo"]
    header_data = header.split(",")
    for type in types:
        out_file = open(os.path.join(out_dir, "%s.txt" % type), "w")
        out_file.write("samples\tcols\tdist\tdata\ttype\n")
        path_idx = get_indices_containing(header_data, type)
        for result in results:
            for idx in path_idx:
                out_file.write("%s\t%s\t%s\t%s\t%s\n" % (result[0], result[1], result[2], result[idx], header_data[idx]))


def get_indices_containing(header, s):
    idx = []
    for i, elem in enumerate(header):
        if s in elem:
            idx.append(i)
    return idx


if __name__ == "__main__":
    main()




