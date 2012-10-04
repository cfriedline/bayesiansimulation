__author__ = 'chris'
import dendropy
import os
import multiprocessing as mp

def chunks(list, size):
    return [list[i:i+size] for i in xrange(0, len(list), size)]

def get_files(dir):
    nex_files = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            f = os.path.join(root, file)
            if f.endswith(".nex"):
                nex_files.append(f)
    return nex_files

def convert_to_phylip(file):
    print "processing", file
    data = dendropy.CharacterMatrix.get_from_path(file, "nexus")
    assert isinstance(data, dendropy.StandardCharacterMatrix)

    out_file = os.path.join(os.path.dirname(file), os.path.basename(file) + ".phy")
    for taxa in data:
        assert isinstance(taxa, dendropy.Taxon)
        taxa.label = taxa.label.ljust(10)
    data.write(open(out_file, "w"), "phylip")


def main():
    files = get_files("/Users/chris/simulations")
    pool = mp.Pool()
    for file in files:
        pool.apply_async(convert_to_phylip, (file,))
    pool.close()
    pool.join()


if __name__ == '__main__':
    main()
    print "Done!"


