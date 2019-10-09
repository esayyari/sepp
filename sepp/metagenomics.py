import os
import tempfile
import re
import sepp
from sepp.alignment import MutableAlignment
from sepp.alignment import _write_fasta
from sepp.config import options
'''
Collection of functions for metagenomic pipeline for taxonomic classification
Created on June 3, 2014

@author: namphuon
'''
global character_map, taxon_map, level_map, key_map, marker_genes, cog_genes
character_map = {'A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'T': 'A',
                 't': 'a', 'G': 'C', 'g': 'c', '-': '-'}
global levels
levels = ["species", "genus", "family", "order", "class", "phylum"]
#marker_genes = [
#    "nusA", "rplB", "rplK", "rplS", "rpsE", "rpsS",  "pgk", "rplC", "rplL",
#    "rplT", "rpsI", "smpB", "dnaG", "pyrg",  "rplD", "rplM", "rpmA", "rpsJ",
#    "frr", "pyrG1", "rplE", "rplN",  "rpsB", "rpsK", "infC", "rplA", "rplF",
#    "rplP", "rpsC", "rpsM"]
marker_genes = ["p0000", "p0004", "p0006", "p0007", "p0009", "p0012", "p0013", "p0014", "p0015", "p0016", "p0018", "p0019", "p0020", "p0022", "p0024", "p0025", "p0026", "p0029", "p0031", "p0032", "p0033", "p0034", "p0035", "p0036", "p0037", "p0038", "p0039", "p0040", "p0041", "p0043", "p0044", "p0046", "p0048", "p0049", "p0050", "p0051", "p0052", "p0054", "p0055", "p0058", "p0059", "p0060", "p0062", "p0063", "p0064", "p0065", "p0066", "p0067", "p0068", "p0069", "p0070", "p0072", "p0073", "p0074", "p0076", "p0077", "p0078", "p0079", "p0080", "p0081", "p0082", "p0083", "p0084", "p0086", "p0087", "p0088", "p0089", "p0090", "p0091", "p0093", "p0094", "p0095", "p0096", "p0097", "p0098", "p0099", "p0101", "p0102", "p0103", "p0104", "p0105", "p0106", "p0107", "p0108", "p0110", "p0112", "p0114", "p0115", "p0116", "p0117", "p0118", "p0119", "p0120", "p0121", "p0122", "p0123", "p0125", "p0126", "p0128", "p0129", "p0130", "p0131", "p0133", "p0134", "p0135", "p0136", "p0137", "p0138", "p0139", "p0140", "p0141", "p0142", "p0143", "p0144", "p0145", "p0146", "p0147", "p0148", "p0149", "p0150", "p0151", "p0152", "p0153", "p0154", "p0155", "p0156", "p0157", "p0159", "p0161", "p0162", "p0163", "p0164", "p0165", "p0166", "p0167", "p0168", "p0169", "p0170", "p0171", "p0172", "p0173", "p0174", "p0175", "p0176", "p0177", "p0178", "p0179", "p0180", "p0181", "p0182", "p0183", "p0184", "p0185", "p0186", "p0187", "p0188", "p0189", "p0190", "p0191", "p0192", "p0193", "p0194", "p0195", "p0196", "p0197", "p0198", "p0199", "p0200", "p0201", "p0202", "p0203", "p0204", "p0205", "p0206", "p0207", "p0208", "p0209", "p0210", "p0211", "p0213", "p0214", "p0215", "p0216", "p0217", "p0219", "p0220", "p0221", "p0222", "p0223", "p0224", "p0225", "p0226", "p0227", "p0228", "p0229", "p0230", "p0231", "p0232", "p0233", "p0234", "p0235", "p0236", "p0237", "p0238", "p0239", "p0241", "p0242", "p0244", "p0246", "p0247", "p0248", "p0249", "p0250", "p0252", "p0253", "p0254", "p0255", "p0256", "p0257", "p0258", "p0259", "p0260", "p0261", "p0262", "p0263", "p0264", "p0265", "p0266", "p0267", "p0268", "p0269", "p0270", "p0271", "p0272", "p0273", "p0274", "p0277", "p0278", "p0279", "p0280", "p0282", "p0283", "p0284", "p0285", "p0286", "p0287", "p0288", "p0289", "p0290", "p0291", "p0292", "p0293", "p0294", "p0295", "p0297", "p0298", "p0299", "p0300", "p0301", "p0303", "p0304", "p0305", "p0306", "p0307", "p0308", "p0309", "p0310", "p0311", "p0312", "p0313", "p0314", "p0315", "p0316", "p0317", "p0318", "p0319", "p0320", "p0321", "p0322", "p0323", "p0324", "p0325", "p0326", "p0328", "p0330", "p0332", "p0333", "p0334", "p0335", "p0337", "p0338", "p0339", "p0340", "p0341", "p0342", "p0343", "p0345", "p0346", "p0347", "p0348", "p0350", "p0351", "p0352", "p0354", "p0355", "p0356", "p0357", "p0358", "p0360", "p0361", "p0362", "p0363", "p0364", "p0365", "p0366", "p0367", "p0368", "p0369", "p0370", "p0371", "p0372", "p0373", "p0375", "p0376", "p0377", "p0378", "p0379", "p0380", "p0381", "p0382", "p0383", "p0384", "p0385", "p0386", "p0387", "p0388", "p0390", "p0391", "p0392", "p0393", "p0394", "p0395", "p0396", "p0397", "p0398", "p0399"]

cog_genes = [
    "COG0049", "COG0088", "COG0094", "COG0100", "COG0184", "COG0201",
    "COG0522", "COG0012", "COG0052", "COG0090", "COG0096", "COG0102",
    "COG0185", "COG0202", "COG0525", "COG0016", "COG0080", "COG0091",
    "COG0097", "COG0103", "COG0186", "COG0215", "COG0533", "COG0018",
    "COG0081", "COG0092", "COG0098", "COG0124", "COG0197", "COG0256",
    "COG0541", "COG0048", "COG0087", "COG0093", "COG0099", "COG0172",
    "COG0200", "COG0495", "COG0552"]


# TODO Fix parameter passing
# TODO Make taxonomy loading a class
def load_taxonomy(taxonomy_file, lower=True):
    global taxon_map, level_map, key_map
    f = open(taxonomy_file, 'r')

    # First line is the keywords for the taxonomy, need to map the keyword to
    # the positional index of each keyword
    results = f.readline().lower().replace('"', '').strip().split(',')
    key_map = dict([(results[i], i) for i in range(0, len(results))])

    # Now fill up taxonomy, level maps keep track of what taxa exist at each
    # level, taxon_map keep track of entire taxonomy
    taxon_map = {}
    level_map = {"species": {}, "genus": {}, "family": {}, "order": {},
                 "class": {}, "phylum": {}}

    for line in f:
        results = line.replace('"', '').strip()
        if (lower):
            results.lower()
        results = results.split(',')
        # insert into taxon map
        taxon_map[results[0]] = results

        # insert into level map
        for level in levels:
            if (results[key_map[level]] == ''):
                continue
            else:
                if (results[key_map[level]] not in level_map[level]):
                    level_map[level][results[key_map[level]]] = {}
                level_map[level][results[key_map[level]]][results[0]] = \
                    results[0]
    return (taxon_map, level_map, key_map)


def build_profile(input, output_directory):
    global taxon_map, level_map, key_map, levels
    temp_dir = tempfile.mkdtemp(dir=options().__getattribute__('tempdir'))
    if (options().bin == 'blast'):
        binned_fragments = blast_to_markers(input, temp_dir)
    else:
        binned_fragments = hmmer_to_markers(input, temp_dir)

    if binned_fragments:
        print("Finished binning")
    else:
        print("Unable to bin any fragments!\n")
        return

    # load up taxonomy for 30 marker genes
    if (options().genes == 'markers'):
        (taxon_map, level_map, key_map) = load_taxonomy(os.path.join(
            options().reference.path, 'refpkg/' + marker_genes[0] + '.refpkg/all_taxon.taxonomy'))
    else:
        (taxon_map, level_map, key_map) = load_taxonomy(os.path.join(
            options().reference.path,
            'refpkg/COG0012.refpkg/all_taxon.taxonomy'))

    # all classifications stored here
    classifications = {}
    classification_files = []
    # Now run TIPP on each fragment
    gene_name = 'sate'
    if (options().genes == 'cogs'):
        gene_name = 'pasta'
    for (gene, frags) in binned_fragments.items():
        # Get size of each marker
        total_taxa = 0
        with open(os.path.join(options().__getattribute__('reference').path,
                  'refpkg/%s.refpkg/%s.size' % (gene, gene_name)), 'r') as f:
            total_taxa = int(f.readline().strip())
        decomp_size = options().alignment_size
        if (decomp_size > total_taxa):
            decomp_size = int(total_taxa / 10)
        cpus = options().cpu
        if (len(frags) < cpus):
            cpus = len(frags)
        extra = ''
        if options().dist is True:
            extra = '-D'
        if options().max_chunk_size is not None:
            extra = extra + '-F %d' % options().max_chunk_size
        if options().cutoff != 0:
            extra = extra+" -C %f" % options().cutoff
        print(
            ('Cmd:\nrun_tipp.py -c %s --cpu %s -m %s -f %s -t %s -adt %s -a '
             '%s -r %s -tx %s -txm %s -at %0.2f -pt %0.2f -A %d -P %d -p %s '
             '-o %s -d %s %s') %
            (options().config_file.name,
             cpus,
             options().molecule,
             temp_dir+"/%s.frags.fas.fixed" % gene,
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/%s.taxonomy' % (gene, gene_name)),
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/%s.tree' % (gene, gene_name)),
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/%s.fasta' % (gene, gene_name)),
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/%s.taxonomy.RAxML_info' % (
                            gene, gene_name)),
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/all_taxon.taxonomy' % gene),
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/species.mapping' % gene),
             options().alignment_threshold,
             0,
             decomp_size,
             total_taxa,
             temp_dir+"/temp_file",
             "tipp_%s" % gene,
             output_directory+"/markers/",
             extra))

        os.system(
            ('run_tipp.py -c %s --cpu %s -m %s -f %s -t %s -adt %s -a %s -r %s'
             ' -tx %s -txm %s -at %0.2f -pt %0.2f -A %d -P %d -p %s -o %s -d '
             '%s %s') %
            (options().config_file.name,
             cpus,
             options().molecule,
             temp_dir+"/%s.frags.fas.fixed" % gene,
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/%s.taxonomy' % (gene, gene_name)),
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/%s.tree' % (gene, gene_name)),
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/%s.fasta' % (gene, gene_name)),
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/%s.taxonomy.RAxML_info' % (
                            gene, gene_name)),
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/all_taxon.taxonomy' % gene),
             os.path.join(options().__getattribute__('reference').path,
                          'refpkg/%s.refpkg/species.mapping' % gene),
             options().alignment_threshold,
             0,
             decomp_size,
             total_taxa,
             temp_dir+"/temp_file",
             "tipp_%s" % gene,
             output_directory+"/markers/",
             extra))
        if (not os.path.exists(output_directory +
                               "/markers/tipp_%s_classification.txt" % gene)):
            continue

        gene_classification = generate_classification(
            output_directory + "/markers/tipp_%s_classification.txt" % gene,
            options().placement_threshold)
        classification_files.append(
            output_directory + "/markers/tipp_%s_classification.txt" % gene)
        # Now write individual classification and also pool classifications
        write_classification(
            gene_classification,
            output_directory + "/markers/tipp_%s.classification" % gene)
        classifications.update(gene_classification)
    remove_unclassified_level(classifications)
    write_classification(classifications,
                         output_directory+"/markers/all.classification")
    write_abundance(classifications, output_directory)

    if (options().dist is True):
        distribution(classification_files, output_directory)


def distribution(classification_files, output_dir):
    global taxon_map, level_map, key_map, levels, level_names
    distribution = {"species": {}, "genus": {}, "family": {}, "order": {},
                    "class": {}, "phylum": {}}
    total_frags = 0
    for class_input in classification_files:
        class_in = open(class_input, 'r')
        frag_info = {"species": {'unclassified': 1},
                     "genus": {'unclassified': 1},
                     "family": {'unclassified': 1},
                     "order": {'unclassified': 1},
                     "class": {'unclassified': 1},
                     "phylum": {'unclassified': 1}}
        old_name = ""
        for line in class_in:
            results = line.strip().split(',')
            if (len(results) > 5):
                results = [results[0], results[1], results[2],
                           results[-2], results[-1]]
            (name, id, rank, probability) = (
                results[0], results[1], results[3], float(results[4]))
            if (rank not in distribution):
                continue
            if (old_name == ""):
                old_name = name
            if (name != old_name):
                total_frags += 1
                assert frag_info['phylum']['unclassified'] != 1
                for clade, cladeval in frag_info.items():
                    for clade_name, cnc in cladeval.items():
                        if (clade_name, cnc not in distribution[clade]):
                            distribution[clade][clade_name] = 0
                        distribution[clade][clade_name] += cnc
                frag_info = {"species": {'unclassified': 1},
                             "genus": {'unclassified': 1},
                             "family": {'unclassified': 1},
                             "order": {'unclassified': 1},
                             "class": {'unclassified': 1},
                             "phylum": {'unclassified': 1}}
                old_name = name
            if (id not in frag_info[rank]):
                frag_info[rank][id] = 0
            frag_info[rank][id] += probability
            frag_info[rank]['unclassified'] -= probability
        total_frags += 1
        assert frag_info['phylum']['unclassified'] != 1
        for clade, cladeval in frag_info.items():
            for clade_name, cnc in cladeval.items():
                if (clade_name not in distribution[clade]):
                    distribution[clade][clade_name] = 0
                distribution[clade][clade_name] += cnc

    level_names = {1: 'species', 2: 'genus', 3: 'family', 4: 'order',
                   5: 'class', 6: 'phylum'}
    for level in level_names:
        f = open(output_dir + "/abundance.distribution.%s.csv" %
                 level_names[level], 'w')
        f.write('taxa\tabundance\n')
        lines = []
        for clade, value in distribution[level_names[level]].items():
            name = clade
            if (name != 'unclassified'):
                name = taxon_map[clade][key_map['tax_name']]
            lines.append('%s\t%0.4f\n' % (name, float(value) / total_frags))
        lines.sort()
        f.write(''.join(lines))
        f.close()
    return distribution


def remove_unclassified_level(classifications, level=6):
    global taxon_map, level_map, key_map, levels
    frags = list(classifications.keys())
    for frag in frags:
        if classifications[frag][level] == 'NA':
            del classifications[frag]


def write_classification(class_input, output):
    '''Writes a classification file
    '''
    class_out = open(output, 'w')
    class_out.write("fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum\n")
    keys = list(class_input.keys())
    keys.sort()
    for frag in keys:
        class_out.write("%s\n" % "\t".join(class_input[frag]))
    class_out.close()


# Fix problem with NA being unclassified
def write_abundance(classifications, output_dir, labels=True,
                    remove_unclassified=True):
    global taxon_map, level_map, key_map, levels

    level_abundance = {
        1: {'total': 0}, 2: {'total': 0}, 3: {'total': 0}, 4: {'total': 0},
        5: {'total': 0}, 6: {'total': 0}}

    level_names = {1: 'species', 2: 'genus', 3: 'family', 4: 'order',
                   5: 'class', 6: 'phylum'}
    for lineage in classifications.values():
        # insert into level map
        for level in range(1, 7):
            if (lineage[level] == 'NA'):
                if ('unclassified' not in level_abundance[level]):
                    level_abundance[level]['unclassified'] = 0
                level_abundance[level]['unclassified'] += 1
                level_abundance[level]['total'] += 1
                # continue
            else:
                if (lineage[level] not in level_abundance[level]):
                    level_abundance[level][lineage[level]] = 0
                level_abundance[level][lineage[level]] += 1
                level_abundance[level]['total'] += 1
    for level in level_names:
        f = open(output_dir + "/abundance.%s.csv" % level_names[level], 'w')
        f.write('taxa\tabundance\n')
        lines = []
        for clade in level_abundance[level]:
            if clade == 'total':
                continue
            name = clade
            if labels and name != 'unclassified':
                name = taxon_map[clade][key_map['tax_name']]
            lines.append('%s\t%0.4f\n' % (
                name, float(level_abundance[level][clade]) / level_abundance[
                    level]['total']))
        lines.sort()
        f.write(''.join(lines))
        f.close()


def generate_classification(class_input, threshold):
    global taxon_map, level_map, key_map, levels
    class_in = open(class_input, 'r')
    level_map_hierarchy = {"species": 0, "genus": 1, "family": 2, "order": 3,
                           "class": 4, "phylum": 5, "root": 6}
    # Need to keep track of last line so we can determine when we switch to
    # new classification
    old_name = ""
    old_probability = 1
    old_id = ""
    old_rank = ""

    # keep track of all fragment names
    names = {}
    classification = {}
    for line in class_in:
        results = line.strip().split(',')
        if (len(results) > 5):
            results = [results[0], results[1], results[2],
                       results[-2], results[-1]]
        (name, id, rank, probability) = (
            results[0], results[1], results[3], float(results[4]))
        names[name] = name
        if (name != old_name):
            # when we switch to new fragment, output last classification for
            # old fragment
            if (old_name != ""):
                lineage = taxon_map[old_id]
                output_line = [old_name]
                for level in levels:
                    clade = lineage[key_map[level]]
                    if (clade == ""):
                        clade = "NA"
                    output_line.append(clade)
                classification[old_name] = output_line
            old_name = name
            old_rank = "root"
            old_probability = 1
            old_id = '1'

        # Switch to new rank if the new probability is higher than threshold
        # and our rank is more specific than our original rank
        if (
            rank in level_map_hierarchy and
            (level_map_hierarchy[old_rank] > level_map_hierarchy[rank]) and
            (probability > threshold)
           ):
            old_rank = rank
            old_probability = probability
            old_id = id
        # Switch to new rank if the new rank matches old rank but has higher
        # probability
        elif (
              rank in level_map_hierarchy and
              (level_map_hierarchy[old_rank] == level_map_hierarchy[rank]) and
              (probability > old_probability)
             ):
            old_rank = rank
            old_probability = probability
            old_id = id

    if old_id in taxon_map:
        lineage = taxon_map[old_id]
        output_line = [old_name]
        for level in levels:
            clade = lineage[key_map[level]]
            if (clade == ""):
                clade = "NA"
            output_line.append(clade)
        classification[name] = output_line
    return classification


def hmmer_to_markers(input, temp_dir):
    global marker_genes
    fragments = MutableAlignment()
    fragments.read_filepath(input)

    reverse = dict([(name+'_rev', reverse_sequence(seq))
                    for (name, seq) in fragments.items()])
    all_frags = MutableAlignment()
    all_frags.set_alignment(fragments)
    all_frags.set_alignment(reverse)
    frag_file = temp_dir+"/frags.fas"
    _write_fasta(all_frags, frag_file)

    # Now bin the fragments
    frag_scores = dict([(name, [-10000, 'NA', 'NA'])
                        for name in fragments.keys()])
    gene_set = marker_genes
    align_name = 'sate'
    if (options().genes == 'cogs'):
        gene_set = cog_genes
        align_name = 'pasta'
    for gene in gene_set:
        # Now run HMMER search
        hmmer_search(
            frag_file,
            os.path.join(
                options().__getattribute__('reference').path,
                'refpkg/%s.refpkg/%.profile' % (gene, align_name)),
            temp_dir + "/%s.out" % gene)
        results = read_hmmsearch_results(temp_dir + "/%s.out" % gene)

        # Now select best direction for each frag
        for name, value in results.items():
            bitscore = value[1]
            direction = 'forward'
            true_name = name
            if (name.find('_rev') != -1):
                true_name = true_name.replace('_rev', '')
                direction = 'reverse'
            if frag_scores[true_name][0] < bitscore:
                frag_scores[true_name] = [bitscore, gene, direction]

    # Now bin the fragments
    genes = dict([])
    for name, val in frag_scores.items():
        if (val[1] not in genes):
            genes[val[1]] = {}
        if (val[2] == 'forward'):
            genes[val[1]][name] = fragments[name]
        else:
            genes[val[1]][name] = reverse_sequence(fragments[name])
    genes.pop("NA", None)
    for gene, seq in genes.items():
        gene_file = temp_dir + "/%s.frags.fas" % gene
        _write_fasta(seq, gene_file + ".fixed")
    return genes


def blast_to_markers(input, temp_dir):
    fragments = MutableAlignment()
    fragments.read_filepath(input)

    if (options().gene is None):
        # First blast sequences against all markers
        blast_results = temp_dir + "/blast.out"
        if (options().blast_file is None):
            print("Blasting fragments against marker dataset\n")
            blast_fragments(input, blast_results)
        else:
            blast_results = options().blast_file
        # Next bin the blast hits to the best gene
        gene_binning = bin_blast_results(blast_results)
    else:
        gene_binning = {options().gene: list(fragments.keys())}
    # Now figure out direction of fragments
    binned_fragments = dict([
        (gene, dict([(seq_name, fragments[seq_name])
                     for seq_name in gene_binning[gene]]))
        for gene in gene_binning])
    print("Finding best orientation of reads\n")
    align_name = 'sate'
    if (options().genes == 'cogs'):
        align_name = 'pasta'
    for (gene, frags) in binned_fragments.items():
        # Add reverse complement sequence
        frags_rev = dict([(name + '_rev', reverse_sequence(seq))
                          for (name, seq) in frags.items()])
        gene_frags = MutableAlignment()
        gene_frags.set_alignment(frags)
        gene_frags.set_alignment(frags_rev)
        gene_file = temp_dir + "/%s.frags.fas" % gene
        _write_fasta(gene_frags, gene_file)

        # Now run HMMER search
        hmmer_search(
            gene_file,
            os.path.join(
                options().__getattribute__('reference').path,
                'refpkg/%s.refpkg/%s.hmm' % (gene, align_name)),
            temp_dir + "/%s.out" % gene)
        results = read_hmmsearch_results(temp_dir + "/%s.out" % gene)

        # Now select best direction for each frag
        for key in frags:
            forward_score = -10000
            backward_score = -10000
            if (key in results):
                forward_score = results[key][1]
            if (key+"_rev" in results):
                backward_score = results[key + "_rev"][1]
            if (backward_score > forward_score):
                frags[key] = gene_frags[key + "_rev"]

        # Now write to file
        _write_fasta(frags, gene_file + ".fixed")
        binned_fragments[gene] = frags
    return binned_fragments


def read_hmmsearch_results(input):
    # Group 1 (e-value) 2 (bitscore) and 9 (taxon name) contain the
    # relevant information, other ones can be ignored unless we plan to do
    # something later
    pattern = re.compile(
        r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)"
        r"\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)")
    start_reading = False
    infile = open(input)
    results = {}
    for line in infile:
        line = line.strip()
        if (not start_reading and line.startswith("E-value") is True):
            start_reading = True
        elif (start_reading and line == ""):
            start_reading = False
            break
        elif (start_reading):
            matches = pattern.search(line)
            if (matches is not None and matches.group(0).find("--") == -1):
                results[matches.group(9).strip()] = (
                    float(matches.group(1).strip()),
                    float(matches.group(2).strip()))
    return results


def read_mapping(input, header=False, delimiter='\t'):
    '''Read a mapping file
    '''
    d = {}
    with open(input) as f:
        for line in f:
            if (header is True):
                next
            results = line.strip().split(delimiter)
            d[results[0]] = results
    return d


def bin_blast_results(input):
    # Map the blast results to the markers
    gene_mapping = read_mapping(
        os.path.join(
            options().__getattribute__('reference').path,
            'blast/%s/seq2marker.tab' % options().genes))

    genes = {}
    with open(input) as f:
        for line in f:
            results = line.split('\t')
            gene = gene_mapping[results[1]][1]
            if gene in genes:
                genes[gene].append(results[0])
            else:
                genes[gene] = [results[0]]
    return genes


def hmmer_search(input, hmmer, output):
    '''Blast the fragments against all marker genes+16S sequences, return
    output'''
    os.system('%s --noali -E 10000 --cpu %d -o %s %s %s' % (
        options().__getattribute__('hmmsearch').path,
        options().cpu, output, hmmer, input))


def blast_fragments(input, output):
    '''Blast the fragments against all marker genes+16S sequences, return
    output'''
    os.system(
        ('%s -db %s -outfmt 6 -query %s -out %s -num_threads %d '
         '-max_target_seqs 1 ') %
        (options().__getattribute__('blast').path,
         os.path.join(
            options().__getattribute__('reference').path,
            "blast/%s/alignment.fasta.db" % options().genes),
         input, output, options().cpu))


def reverse_sequence(sequence):
    global character_map
    #  Reverse complement the sequence
    return "".join([character_map.get(a, a) for a in sequence[::-1]])


def augment_parser():
    # default_settings['DEF_P'] = (100 , "Number of taxa (i.e.
    # no decomposition)")
    parser = sepp.config.get_parser()

    tippGroup = parser.add_argument_group(
        "TIPP Options".upper(),
        "These arguments set settings specific to TIPP")

    tippGroup.add_argument(
        "-at", "--alignmentThreshold", type=float,
        dest="alignment_threshold", metavar="N",
        default=0.0,
        help="Enough alignment subsets are selected to reach a commulative "
             "probability of N. "
             "This should be a number between 0 and 1 [default: 0.95]")

    tippGroup.add_argument(
        "-pt", "--placementThreshold", type=float,
        dest="placement_threshold", metavar="N",
        default=0.0,
        help="Enough placements are selected to reach a commulative "
             "probability of N. "
             "This should be a number between 0 and 1 [default: 0.95]")
    tippGroup.add_argument(
        "-g", "--gene", type=str,
        dest="gene", metavar="N",
        default=None,
        help="Classify on only the specified gene. ")

    tippGroup.add_argument(
        "-b", "--blast_file", type=str,
        dest="blast_file", metavar="N",
        default=None,
        help="Blast file with fragments already binned. ")

    tippGroup.add_argument(
        "-bin", "--bin_using", type=str,
        dest="bin", metavar="N",
        default="blast",
        help="Tool for binning")

    tippGroup.add_argument(
        "-D", "--dist",
        dest="dist", action='store_true',
        default=False,
        help="Treat fragments as distribution")

    tippGroup.add_argument(
        "-C", "--cutoff", type=float,
        dest="cutoff", metavar="N",
        default=0.0,
        help="Placement probability requirement to count toward the "
             "distribution. "
             "This should be a number between 0 and 1 [default: 0.0]")

    tippGroup.add_argument(
        "-G", "--genes", type=str,
        dest="genes", metavar="GENES",
        default='markers',
        help="Use markers or cogs genes [default: markers]")


def main():
    augment_parser()
    sepp.config._options_singelton = sepp.config._parse_options()
    if (options().alignment_size is None):
        options().alignment_size = int(total_taxa / 10)
    input = options().fragment_file.name
    output_directory = options().outdir
    build_profile(input, output_directory)


if __name__ == '__main__':
    main()
