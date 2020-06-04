"""
Convert a pedigree "fam" file into a "def" file for input into Ped-sim.

The input file can be tab or space separated and should have a header with at
least the following columns in this order:

ID FATHER MOTHER SEX

The IDs can be either integers or strings, and the sex should be 1 for male, 2
for female.

Author: Sara Mathieson
Date: 05/08/20
"""

from collections import defaultdict
import optparse
import sys

################################################################################
# CLASSES
################################################################################

class Generation:
    """Keeps track of the ordered branches in a generation"""

    def __init__(self):
        self.branches = []

    def __str__(self):
        return ", ".join(self.branches)

    def add_branch(self, bid):
        """Add a branch with the given individual's ID"""
        assert bid not in self.branches
        self.branches.append(bid)

    def find_branch(self, pid):
        """Given parent ID, find the branch index of the parent"""
        lst = []
        for i, bid in enumerate(self.branches):
            if pid == bid:
                lst.append(i)

        # make sure parent appears in at most one branch
        assert len(lst) <= 1
        if len(lst) == 1:
            return lst[0]
        return None

class Individual:
    """Keeps track information about each individual in the pedigree"""

    def __init__(self, iid, father, mother, sex):
        self.iid = iid
        self.spouses = set()
        self.father = father
        self.mother = mother
        self.sex = sex
        self.children = []
        self.founder = False

    def add_spouse(self, sid):
        """Add a spouse (note the type is a set)"""
        self.spouses.add(sid)

    def add_child(self, cid):
        """Add a child"""
        assert cid not in self.children
        self.children.append(cid)

################################################################################
# PARSE ARGS and INPUT FILE
################################################################################

def parse_args():
    """Parse command line arguments."""
    parser = optparse.OptionParser(description='fam to def file for Ped-sim')

    parser.add_option('-i', '--ped_filename', type='string', \
        help='path to input fam file')
    parser.add_option('-o', '--def_filename', type='string', \
        help='path to output def file')
    parser.add_option('-n', '--name', type='string', default='ped', \
        help='name of the pedigree to use in the def file')

    (opts, args) = parser.parse_args()

    mandatories = ['ped_filename', 'def_filename',]
    for m in mandatories:
        if not opts.__dict__[m]:
            print('mandatory option ' + m + ' is missing\n')
            parser.print_help()
            sys.exit()

    return opts

def read_indvs(filename):
    """Parse input fam file to create Individuals"""
    all_indvs = {} # key=ID (string), val=Individual
    num_founder = 0

    # read file
    f = open(filename,'r')
    f.readline() # header
    for line in f:
        tokens = line.strip().split()
        indv = Individual(tokens[0], tokens[1], tokens[2], tokens[3])
        all_indvs[tokens[0]] = indv
    f.close()

    # process all the individuals to add spouses and children
    for iid, indv in all_indvs.items():
        fid = indv.father
        mid = indv.mother

        # non-founder individual
        if fid != '0' and mid != '0':
            all_indvs[fid].add_spouse(mid)
            all_indvs[mid].add_spouse(fid)
            all_indvs[fid].add_child(iid)
            all_indvs[mid].add_child(iid)
        else:
            indv.founder = True
            num_founder += 1

    print("fam file, num individuals:", len(all_indvs))
    print("fam file, num founders:", num_founder)
    return all_indvs

################################################################################
# CREATE GENERATIONS
################################################################################

def all_generations(all_indvs):
    """From dictionary of all individuals, iteratively create generations"""
    all_gens = []
    used_ids = [] # keep track of which individuals we've processed

    # set up first generation with founders
    prev_gen = Generation()
    for iid, indv in all_indvs.items():
        if indv.founder:

            # one spouse
            if len(indv.spouses) == 1:
                sid = list(indv.spouses)[0]
                if all_indvs[sid].founder:
                    if indv.sex == "1":
                        assert iid not in used_ids and sid not in used_ids
                        prev_gen.add_branch(iid)
                        used_ids.append(iid)
                        used_ids.append(sid)

            # multiple spouses
            else:
                assert iid not in used_ids
                prev_gen.add_branch(iid)
                used_ids.append(iid)

    # continue creating generations while we still have active branches
    while len(prev_gen.branches) > 0:
        all_gens.append(prev_gen)
        next_gen = Generation()

        # find all the children of all parents in the previous generation
        for bid in prev_gen.branches:
            parent = all_indvs[bid]

            for cid in parent.children:
                child = all_indvs[cid]
                fid = child.father
                mid = child.mother
                father = all_indvs[fid]
                mother = all_indvs[mid]

                # neither parent is founder
                if not father.founder and not mother.founder:
                    # make sure parents are NOT in *current* generation
                    if fid in used_ids and mid in used_ids and fid not in \
                        next_gen.branches and mid not in next_gen.branches:

                        # add a branch for this child
                        if cid not in used_ids:
                            next_gen.add_branch(cid)
                            used_ids.append(cid)

                    # WAIT to add until generation after one of the parents
                    else:
                        pass

                # one parent is founder
                else:
                    # pass over married-in with multiple spouses
                    if parent.founder and len(parent.spouses) > 1:
                        pass
                    else:
                        if fid == parent.iid:
                            assert mother.founder
                            if mid not in used_ids:
                                used_ids.append(mid)
                        else:
                            assert father.founder
                            if fid not in used_ids:
                                used_ids.append(fid)

                        # add a branch for this child
                        if cid not in used_ids:
                            next_gen.add_branch(cid)
                            used_ids.append(cid)

        # set up for the next generation
        prev_gen = next_gen

    used_ids = set(used_ids)
    print("def file, num individuals:", len(used_ids))

    # look at number of individuals in branches
    in_branches = []
    for gen in all_gens:
        in_branches.extend(gen.branches)
    in_branches = set(in_branches)

    # founders are counted in branches in gen 0
    num_founders = len(all_indvs) - len(in_branches) + len(all_gens[0].branches)
    return all_gens

def find_gen_branch(iid, g, all_gens):
    """
    Based on the ID of an individual and the current generation, search through
    previous generations for the branch index.
    """

    for prev_gen in range(g-1, -1, -1):
        b = all_gens[prev_gen].find_branch(iid) # indexed from zero
        if b != None:
            return prev_gen, b # return as soon as we find it

    return None, None

def pretty_print(s_dict, all_indvs):
    """From a dictionary of children for each couple, create the "def" format"""
    s = " "
    founder_ids = []

    # go through each couple
    for k, v in s_dict.items():

        branches = []
        string = None
        # all children should have same parents here
        for item in v:
            branches.append(item[0])
            if string == None:
                string = item[1]
            else:
                assert item[1] == string

        fid = k.split('+')[0]
        mid = k.split('+')[1]
        if all_indvs[fid].founder:
            if fid not in founder_ids:
                founder_ids.append(fid)
        if all_indvs[mid].founder:
            if mid not in founder_ids:
                founder_ids.append(mid)

        # one or both spouses is a founder (this is just sanity checking)
        if "_" not in string:
            assert all_indvs[fid].founder or all_indvs[mid].founder
        # check for married-in with multiple spouses
        elif all_indvs[fid].founder or all_indvs[mid].founder:
            #print(fid, mid, 'married-in with multiple spouses')
            pass

        # sort children in each branch
        branches.sort()

        # one child
        if len(branches) == 1:
            s += str(branches[0]) + ":" + string + " "
        # children are consecutive
        elif branches[-1] - branches[0] + 1 == len(branches):
            s += str(branches[0]) + "-" + str(branches[-1]) + ":" + string + " "
        # children are not consecutive
        else:
            s += ",".join([str(branches[i]) for i in range(len(branches))]) + \
                ":" + string + " "

    return s, founder_ids

################################################################################
# MAIN
################################################################################

def main():
    """Orchestrate generation creation and writing of output file"""
    opts = parse_args()
    all_indvs = read_indvs(opts.ped_filename)
    all_gens = all_generations(all_indvs)
    founder_ids = []

    # set up def file
    def_file = open(opts.def_filename,'w')
    def_file.write("def " + opts.name + " 1 " + str(len(all_gens)) + "\n")

    # write each generation in def format
    for g in range(len(all_gens)):
        gen = all_gens[g]

        # keeps track of children with the same two parents
        s_dict = defaultdict(list)
        if g != 0:

            # go through all the branches to find parents in the prev generation
            for b, bid in enumerate(gen.branches):
                assert not all_indvs[bid].founder

                # IDs of father and mother
                fid = all_indvs[bid].father
                mid = all_indvs[bid].mother

                # NOTE: if founder, omit spouse to create new founder spouse
                # NOTE: one parent must be from the previous generation
                # find (generation, branch_index) of indv
                fg, fb = find_gen_branch(fid, g, all_gens)
                mg, mb = find_gen_branch(mid, g, all_gens)

                # CASE 1: founder father
                if fg == None:
                    assert mg == g-1 # mother in prev generation
                    assert all_indvs[fid].founder
                    add_str = str(mb+1) # just include branch index

                # CASE 2: founder mother
                elif mg == None:
                    assert fg == g-1 # father in prev generation
                    assert all_indvs[mid].founder
                    add_str = str(fb+1) # just include branch index

                # CASE 3: both father and mother in prev generation
                elif fg == g-1 and mg == g-1:
                    add_str = str(fb+1) + "_" + str(mb+1)

                # CASE 4: father in prev generation
                elif fg == g-1:
                    add_str = str(fb+1) + "_" + str(mb+1) + "^" + str(mg+1)

                # CASE 5: mother in prev generation
                elif mg == g-1:
                    add_str = str(mb+1) + "_" + str(fb+1) + "^" + str(fg+1)

                # should not get here!
                else:
                    print("ERROR: no valid parents for indv", bid)

                # add this indv to correct parents
                parent_str = fid+"+"+mid
                s_dict[parent_str].append([b+1, add_str])

        # convert to def format
        s, f_ids_gen = pretty_print(s_dict, all_indvs)

        # keep track of founder IDs for debugging
        founder_ids.extend(f_ids_gen)

        # write to def
        def_file.write(str(g+1) + " 1 " + str(len(gen.branches)) + s + "\n")

    founder_ids = set(founder_ids)
    print("def file, num founder:", len(founder_ids))

if __name__ == "__main__":
    main()
