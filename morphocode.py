__author__ = 'chris'
import app
import dendropy
import numpy
import copy

def read_source_file(file):
    f = open(file, "r")
    line_num = 0
    taxa = []
    data = []
    for line in f:
        line =line.strip()
        if line_num > 0:
            d = line.split("\t")
            taxa.append(d[0])
            data.append(d[1:])
        line_num += 1
    return taxa, data

def read_result_file(file):
    return dendropy.StandardCharacterMatrix.get_from_path(file, "nexus")

def get_column(data, col):
    return [row[col] for row in data]

def get_gap_weight(value, states, min, max):
    if value == '-':
        return value
    else:
        return int(round(((float(value) - min) * (states - 1)) / (max - min)))

def get_min_max(col):
    col = copy.copy(col)
    count = col.count('-')
    for i in xrange(count):
        col.remove('-')
    c = [float(elem) for elem in col]
    return min(c), max(c)



def get_gap_weights(data):
    cols = len(data[0])
    gaps = []
    for i in xrange(cols):
        col = get_column(data, i)
        min_max = get_min_max(col)
        col_gaps = []
        for val in col:
            col_gaps.append(get_gap_weight(val, 10, min_max[0], min_max[1]))
        gaps.append(col_gaps)
    gaps = numpy.array(gaps)
    return gaps.transpose()

def main():
    taxa, data = read_source_file("Sample tab-delimited text.txt")
    gaps = get_gap_weights(data)
    coded = read_result_file("Morphometric data.nex")
    my_gaps = []
    morpho_gaps = []

    for row in gaps:
        my_gaps.append(''.join(row))

    for taxon in coded.taxon_set:
        morpho_gaps.append(coded[taxon])

    for i in xrange(len(my_gaps)):
        assert my_gaps[i] == str(morpho_gaps[i])
        print "%s (me) = %s (morphocode)" % (my_gaps[i], morpho_gaps[i])



if __name__ == '__main__':
    main()
    print "Done!"


