# 
# title:   asteroids.py
# summary: Implementation of Blackledge's ASTEROIDS algorithm
# author: 
#          Yue Zhang (yue.zhang at nih.gov)
# date:    2/20/2021
#
import sys

SZ_WORLD = 1024 * 32  # Total number of structures (1e6)
N_PTS = 64  # Number of "RDC" data points
N_ENSEMBLE = 128  # Number of ensembles per. generation (100)


# Note for OS: setting can be changed 
N_GENERATIONS = 200  # Number of generations (2000)

RANDOM_MUTATION_RATE = 0.01  # Random mutation rate (1%)

MAX_COEFF = 5.0  # "Toy" data: Maximum coefficient
MAX_K = float(N_PTS)  # "Toy" data: Maximum wavenumber

DEBUG = 2  # Debug level

import random, math, sys, os
import numpy as np

# random.seed(1)

###########################################################################
# General functions
###########################################################################

def read_error_weights(fn=''):
    error, weight_factor = [], []
    with open(fn) as f:
        lines = f.readlines()
    for l in lines:
        l = l.split()
        if l:
            error.append(float(l[1]))
            weight_factor.append(float(l[2]))

    return error, weight_factor



def generate_world(path=""):
    global SZ_WORLD, N_PTS

    world = []
    index = []

    # If no path given, generate a random set of sine and cosine functions
    if not path:
        for i in range(SZ_WORLD):
            if i % 1024 == 0:
                dbg("Generating set: %10i of %10i" % (i, SZ_WORLD), 0)

            c = random.uniform(MAX_COEFF, -MAX_COEFF)
            k = random.uniform(0.0, MAX_K)

            if random.uniform(0.0, 1.0) > 0.5:
                func = math.cos
                fstr = "cos"
            else:
                func = math.sin
                fstr = "sin"

            r = []
            for j in range(N_PTS):
                r.append(c * func(2.0 * math.pi * k * float(j) / MAX_K))

            world.append(np.array(r))
            index.append(i)

        return world, index

    npath = os.path.abspath(os.path.expanduser(os.path.expandvars(path)))

    if not os.path.exists(npath):
        sys.stderr.write("Unable to find path: %s\n" % path)
        sys.exit(1)

    # If the path is given, assume it contains data (in one file or many
    # files) corresponding to calculated data for each PDB file; the data
    # is assumed to be two columns; the first (residue number) is ignored
    # here, because it's not needed.  If the world and target functions
    # follow the same order, the program does not need to know that
    # the entry is 1JNH for residue 1, the second is 1JCAHA for residue 1,
    # etc.  The data can be entered in any order provided the order is
    # consistent in the world and target files.

    # First option: a directory is given, containing many data files
    if os.path.isdir(npath):
        # could be modified to recurse, filter, etc.
        files = ["%s/%s" % (npath, fn) for fn in os.listdir(npath)]
        n_files = 0

        # dry run; determine the length of the dataset
        f = open(files[0])
        l = f.readline()
        ndata = 0
        while l:
            l = l.strip()
            if not l or l[0] == "#":
                l = f.readline()
                continue
            ndata = ndata + 1
            l = f.readline()
        f.close()
        N_PTS = ndata
        print(N_PTS)
        for i in range(len(files)):
            if i % 1024 == 0:
                dbg("Reading file: %10i of %10i" % (i, len(files)), 0)

            try:
                fn = files[i]
                f = open(fn)
                l = f.readline()
                data = []

                while l:
                    l = l.strip()
                    if not l or l[0] == "#":
                        l = f.readline()
                        continue
                    toks = l.split()
                    data.append(float(toks[1]))  # ignore 1st column
                    l = f.readline()
                f.close()
            except:
                dbg("! Error reading file: %s" % fn, 0)
                continue
                # raise

            if not len(data) == N_PTS:
                print(len(data), N_PTS)
                dbg("! File doesn't contain %s points: %s" % (N_PTS, fn), 0)
                continue

            world.append(np.array(data))
            index.append(fn)
            n_files = n_files + 1

        SZ_WORLD = n_files
        return world, index

    # Second option: a single file is given, containing all the datasets
    # with each dataset separated by a blank line.

    # Determine the size of each dataset
    f = open(npath)
    l = f.readline()

    ndata = 0

    while l:
        # print l[:-1]
        l = l.strip()
        if not l:
            break
        if l[0] == "#":
            l = f.readline()
            continue
        ndata = ndata + 1
        l = f.readline()
    f.close()
    N_PTS = ndata

    # Parse the file for all the datasets
    f = open(npath)
    l = f.readline()

    data = []
    didx = 0

    while l:
        # print l[:-1]
        l = l.strip()
        if not l:
            if not len(data) == 0:
                if didx % 1024 == 0:
                    dbg("Adding dataset: %10i" % didx, 0)

                assert len(data) == N_PTS
                world.append(np.array(data))
                index.append(didx)
                didx = didx + 1
                data = []

            l = f.readline()
            continue

        if l[0] == "#":
            l = f.readline()

            continue

        toks = l.split()
        data.append(float(toks[1]))  # ignore 1st column
        l = f.readline()

    if not len(data) == 0:
        assert len(data) == N_PTS
        world.append(np.array(data))
        index.append(didx)
        didx = didx + 1
        data = []

    f.close()
    SZ_WORLD = didx

    return world, index


def generate_target(path=""):
    result = []

    if not path:
        for i in range(N_PTS):
            result.append(
                0.5 * math.sin(2.0 * math.pi / 12.0 * float(i))
                + 0.5 * math.cos(2.0 * math.pi / 49.0 * float(i))
            )
        return np.array(result)

    npath = os.path.abspath(os.path.expanduser(os.path.expandvars(path)))

    if not os.path.exists(npath):
        sys.stderr.write("Unable to find target path: %s\n" % path)
        sys.exit(1)

    f = open(npath)
    l = f.readline()
    data = []

    while l:
        l = l.strip()
        if not l or l[0] == "#":
            l = f.readline()
            continue
        toks = l.split()
        data.append(float(toks[1]))  # ignore 1st column
        l = f.readline()
    f.close()

    if not len(data) == N_PTS:
        dbg("! Target doesn't contain %s points (%i)" % (N_PTS, len(data)), 0)
        sys.exit(1)

    return np.array(data)


def factors(n):
    """return the factors of a number"""
    factors = []
    n = int(n)

    for i in range(1, n + 1):
        if not n % i:
            factors.append(i)

    return factors


def schedule_tournaments(length=N_GENERATIONS):
    """create a list of N_GENERATIONS integers, scaling down tournaments"""
    initial = float(N_ENSEMBLE)
    final = float(1)
    delta = float(length)

    possible = factors(N_ENSEMBLE)
    possible.reverse()
    # method of increasing selective pressure is not clear in the
    # ASTEROIDS paper; assuming a linear progression from *initial* to
    # *final*.

    m = (final - initial) / delta
    schedule = []

    for i in range(length):
        linear = max(int(initial + i * m), 1)
        for j in range(1, len(possible)):
            this_fac = possible[j]
            last_fac = possible[j - 1]

            if linear >= this_fac:
                dist_this = abs(this_fac - linear)
                dist_last = abs(last_fac - linear)

                if dist_this < dist_last:
                    schedule.append(this_fac)
                else:
                    schedule.append(last_fac)
                break

    return schedule


def dbg(str, level=0):
    if level <= DEBUG:
        sys.stdout.write("%s\n" % str)
        sys.stdout.flush()


def plot(data, out=sys.stdout):
    close_out = 0
    if type(out) == type(""):
        out = open(out, "w")
        close_out = 1

    w = out.write
    for i in range(len(data)):
        for j in range(len(data[i])):
            w("%5i %13.5e\n" % (j, data[i][j]))

        out.write("\n")

    if close_out:
        out.close()


def safe_dir(a=""):

    if not os.path.exists("%s.dat" % a):
        return "%s.dat" % a
    num = 0
    new_a = "%s_%03i.dat" % (a, num)

    while os.path.exists(new_a):
        num += 1
        new_a = "%s_%03i" % (a, num)
    return new_a


###########################################################################
# Class Definitons
###########################################################################

### Ensemble Class #######################################################


class Ensemble:
    def __init__(self, world, tgt, my_list=None):
        """constructor: make a new Ensemble from a list"""

        self.world = world
        self.target = tgt
        self.structs = []
        self.used = {}
        self.my_score = None

        self.rsize = min(SZ_ENSEMBLE * RANDOM_MUTATION_RATE, 1)

        if my_list:
            assert len(my_list) == SZ_ENSEMBLE
            self.structs = my_list
            for i in my_list:
                assert not self.used.__contains__(i)
                self.used[i] = 1

    def random_unused(self):
        """grab an unused structure for inclusion in the ensemble"""
        struct = int(random.uniform(0, SZ_WORLD))
        n_attempts = 0

        while self.used.__contains__(struct):
            n_attempts = n_attempts + 1
            assert n_attempts < SZ_WORLD
            struct = int(random.uniform(0, SZ_WORLD))

        self.used[struct] = 1
        return struct

    def random(self):
        """generate a new ensemble by random evolution"""
        child = Ensemble(self.world, self.target)

        for i in range(SZ_ENSEMBLE):
            child.structs.append(child.random_unused())
        return child

    def copy(self):
        """copy constructor"""
        new = Ensemble(self.world, self.target)
        new.structs = self.structs[:]
        new.used = self.used.copy()
        new.my_score = self.my_score

        return new

    def mutation(self, internal={}):
        """mutate the current ensemble"""
        n_chose = 0
        gone = {}
        keys = internal.keys()

        child = self.copy()

        while n_chose < child.rsize:
            i = int(random.uniform(0, SZ_ENSEMBLE))
            if gone.__contains__(i):
                continue

            gone[i] = 1
            n_chose = n_chose + 1

        for i in gone.keys():
            del child.used[child.structs[i]]

            # internal mutation
            if len(internal):
                new = random.choice(list(keys))
                n_attempts = 0

                while child.used.__contains__(new):
                    new = random.choice(list(keys))
                    n_attempts = n_attempts + 1

                    if n_attempts > 2 * len(keys):
                        new = child.random_unused()
                        break

                child.used[new] = 1

            # external mutation
            else:
                new = child.random_unused()

            child.structs[i] = new

        child.my_score = None
        return child

    def cross(self, mate):
        """cross with another ensemble, returning an entirely new ensemble"""
        total = {}

        for i in range(len(self.structs)):
            if total.__contains__(self.structs[i]):
                continue
            total[self.structs[i]] = 1

        for i in range(len(mate.structs)):
            if total.__contains__(mate.structs[i]):
                continue
            total[mate.structs[i]] = 1
        # print (len(total), SZ_ENSEMBLE, '$$$$$$$')
        assert len(total) >= SZ_ENSEMBLE

        structs = []
        used = {}

        for i in range(SZ_ENSEMBLE):
            available = total.keys()
            struct = random.choice(list(available))

            structs.append(struct)
            used[struct] = 1
            del total[struct]

        child = Ensemble(self.world, self.target)
        child.structs = structs
        child.used = used

        return child

    def average(self):
        """compute the average of all structures in the ensemble"""
        world = self.world
        return sum([world[i] for i in self.structs]) / SZ_ENSEMBLE

    def score(self):
        """calculate the score compared to the target function"""
        world = self.world
        tgt = self.target

        if self.my_score == None:
            avg = self.average()
            ens = sum([world[i] for i in self.structs]) / SZ_ENSEMBLE
            dif = (avg - tgt) * np.array(WEIGHT_FACTOR)
            dif = dif/tgt
            self.my_score = np.dot(dif, dif)

        return self.my_score

    def plot_all(self, comment=None, out=sys.stdout):
        """plot all structures in an enesmble in a single file"""
        world = self.world

        close_out = 0
        if type(out) == type(""):
            out = open(out, "w")
            close_out = 1

        w = out.write
        if comment:
            w("# %s\n" % comment)

        for i in range(len(self.structs)):
            for j in range(len(world[self.structs[i]])):
                w("%5i %13.5e\n" % (j, world[self.structs[i]][j]))
            w("\n")

        if close_out:
            out.close()
            return
        w("\n")
        out.flush()

    def plot_avg(self, comment=None, out=sys.stdout):
        """plot the average of an ensemble in a single file"""

        close_out = 0
        if type(out) == type(""):
            out = open(out, "w")
            close_out = 1

        avg = self.average()

        w = out.write
        if comment:
            w("# %s\n" % comment)

        for j in range(len(avg)):
            w("%5i %13.5e\n" % (j, avg[j]))

        if close_out:
            out.close()
            return
        w("\n")
        out.flush()

    def cross_reference(self, index, comment=None, out=sys.stdout):
        """print a cross-reference file containing indices"""
        close_out = 0
        if type(out) == type(""):
            out = open(out, "w")
            close_out = 1

        w = out.write
        if comment:
            w("# %s\n" % comment)

        for i in range(len(self.structs)):
            w("%5i %s\n" % (i, index[self.structs[i]]))

        if close_out:
            out.close()
            return
        w("\n")
        out.flush()


### Generation Class ######################################################


class Generation:
    def __init__(self, world, tgt, ensembles=None):
        self.world = world
        self.target = tgt
        self.ensembles = []
        self.rank = []

        if ensembles:
            assert len(ensembles) == N_ENSEMBLE
            self.ensembles = ensembles

            for i in range(N_ENSEMBLE):
                self.rank.append(e.score(), i)

            # self.rank.sort()

    def random(self):
        """create an entire generation by random mutation"""
        new = Generation(self.world, self.target)

        for i in range(N_ENSEMBLE):
            member = Ensemble(new.world, new.target).random()
            new.ensembles.append(member)
            new.rank.append((member.score(), i))

        # new.rank.sort()
        return new

    def mutation(self, internal={}):
        """create a new generation, based on self, by mutation"""
        new = Generation(self.world, self.target)

        for i in range(N_ENSEMBLE):
            old_e = self.ensembles[i]
            new_e = old_e.mutation(internal)
            new.ensembles.append(new_e)
            new.rank.append((new_e.score(), i))

        # new.rank.sort()
        return new

    def copy(self):
        """copy construtor for a generation"""
        new = Generation(self.world, self.target)
        new.ensembles = self.ensembles[:]
        new.rank = self.rank[:]

    def cross(self):
        """create a new generation by crossing"""
        new = Generation(self.world, self.target)

        for i in range(N_ENSEMBLE):
            par1 = int(random.uniform(0, N_ENSEMBLE))
            par2 = par1

            while par2 == par1:
                par2 = int(random.uniform(0, N_ENSEMBLE))

            child = self.ensembles[par1].cross(self.ensembles[par2])
            new.ensembles.append(child)
            new.rank.append((child.score(), i))

        # new.rank.sort()
        return new

    def check_repeats(self):
        """check to see if an ensemble is repeated in a set"""
        table = {}

        for i in range(len(self.ensembles)):
            my_list = self.ensembles[i].structs[:]
            my_list.sort()

            key = tuple(my_list)
            if table.__contains__(key):
                dbg("Ensemble %i is a duplicate with %s" % (i, "table[key]"))
                table[key].append(i)
                return 1
            else:
                table[key] = [i]

        return 0

    def all_structs(self):
        """return a dictionary of all ensembles used in this generation"""
        result = {}

        for i in range(N_ENSEMBLE):
            for j in range(SZ_ENSEMBLE):
                result[self.ensembles[i].structs[j]] = 1

        return result


# def gen_error(target):
#     # error = 5.0
#     error = np.array(EXP_ERRO)

#     x = np.random.sample(len(target))
#     error_factor = ( x * 2.0 * error - error) # random error from -error, error in 0.xxx

#     return target * (1.0 + error_factor)

def gen_error(target, error=15):
    error_factor = (np.random.sample(len(target)) * 2.0 * error - error)/ 100.0 # random error from -error, error in 0.xxx
    return target * (1.0 + error_factor)

###########################################################################
# Main Function
###########################################################################


def main(wpath=None, tpath=None):
    # Generate random data.  This data corresponds to the set of RDCs
    # calculated for many individual files in the entire simulation.
    # Currently SZ_WORLD arrays are generated, each containing N_PTS
    # points.

    world, index = generate_world(wpath)
    target = generate_target(tpath)

    # FOR OS  add for the random error :
    #target = gen_error(target)
     
    out_tgt = safe_dir("target.dat")
    plot([target], out=out_tgt)

    # Create the first generation
    gen = Generation(world, target).random()

    # Determine the number of tournaments in each generation
    schedule = schedule_tournaments()
    dbg("tournament schedule: %s" % "schedule", 4)

    out_avg = safe_dir("best_avg")
    out_ens = safe_dir("best_ens")
    out_idx = safe_dir("best_idx")

    for genidx in range(N_GENERATIONS):
        dbg("\nStarting Generation %i:" % genidx, 2)
        n_tournaments = schedule[genidx]
        n_winners = int(N_ENSEMBLE / n_tournaments)
        total_ensembles = N_ENSEMBLE * 5
        bracket_size = int(total_ensembles / n_tournaments)

        dbg("* n_tournaments   = %i" % n_tournaments, 2)
        dbg("* n_winners       = %i" % n_winners, 2)
        dbg("* total_ensembles = %i" % total_ensembles, 2)
        dbg("* bracket_size    = %i" % bracket_size, 2)

        assert total_ensembles % n_tournaments == 0

        # Debugging purposes only: check for duplicates / inbreeding
        if DEBUG >= 3:
            dbg("- checking duplicates", 3)
            if gen.check_repeats():
                dbg("! parent generation has duplicates!", 3)
                return

        # Generate the next set of ensembles
        rndmut = gen.random()
        intmut = gen.mutation(gen.all_structs())
        extmut = gen.mutation()
        cross = gen.cross()

        if DEBUG >= 50:
            if rndmut.check_repeats():
                dbg("! rndmut generation has duplicates!", 3)
                return
            if intmut.check_repeats():
                dbg("! intmut generation has duplicates!", 3)
                return
            if extmut.check_repeats():
                dbg("! extmut generation has duplicates!", 3)
                return
            if cross.check_repeats():
                dbg("! cross generation has duplicates!", 3)
                return

        # Create a proxy for all players in all the tournaments
        # head to head matches are determined by shuffling this
        # list.
        players = [gen, rndmut, intmut, extmut, cross]
        all_ensembles = []
        for i in range(len(players)):
            for j in range(N_ENSEMBLE):
                score, idx = players[i].rank[j]
                all_ensembles.append((score, i, idx))

        random.shuffle(all_ensembles)

        winners = []
        low_score = None
        low_idx = None
        best_scores = []

        # Complete the tournaments
        for tnum in range(n_tournaments):
            dbg("- running tournment %i" % tnum, 3)

            bracket = all_ensembles[1:bracket_size]
            all_ensembles = all_ensembles[bracket_size:]

            bracket.sort()

            if low_score == None or bracket[0][0] < low_score:
                low_score = bracket[0][0]
                low_idx = len(winners)

            for i in range(n_winners):
                best_scores.append(bracket[i][0])
                winners.append(bracket[i])

        avg_score = sum(best_scores) / len(best_scores)
        dbg("- %i winners selected" % len(winners), 2)
        dbg("- Best Score, round %i: %13.5e" % (genidx, low_score), 2)
        dbg("- Avg. Score, round %i: %13.5e" % (genidx, avg_score), 2)

        # Make a new generation from the previous winners
        next_gen = Generation(world, target)
        enshash = {}

        for i in range(len(winners)):
            score, pidx, eidx = winners[i]
            winsemble = players[pidx].ensembles[eidx]
            slist = winsemble.structs[:]
            slist.sort()
            key = tuple(slist)

            while enshash.__contains__(key):
                pidx0, eidx0, i0 = enshash[key]
                dbg(
                    "- performing external mutation to avoid duplicates\n"
                    "  e1: gen id= %4i, edix= %4i, widx= %4i sc= %13.5e\n"
                    "  e2: gen id= %4i, edix= %4i, widx= %4i sc= %13.5e"
                    % (
                        pidx0,
                        eidx0,
                        i0,
                        players[pidx0].ensembles[eidx0].score(),
                        pidx,
                        eidx,
                        i,
                        players[pidx].ensembles[eidx].score(),
                    ),
                    0,
                )

                winsemble = winsemble.mutation()
                slist = winsemble.structs[:]
                slist.sort()
                key = tuple(slist)

            enshash[key] = (pidx, eidx, i)
            next_gen.ensembles.append(winsemble)
            next_gen.rank.append((winsemble.score(), i))

        gen = next_gen

        best_ensemble = gen.ensembles[low_idx]
        best_ensemble.plot_avg(
            out=out_avg,
            comment="Round %i average (%13.5e)" % (genidx, best_ensemble.score()),
        )
        best_ensemble.plot_all(
            out=out_ens,
            comment="Round %i average (%13.5e)" % (genidx, best_ensemble.score()),
        )
        best_ensemble.cross_reference(
            index,
            out=out_idx,
            comment="Round %i members (%13.5e)" % (genidx, best_ensemble.score()),
        )


if __name__ == "__main__":
    try:
        world = sys.argv[1]
        target = sys.argv[2]
        error_wf  = sys.argv[3]
        SZ_ENSEMBLE = int(sys.argv[-1])  # Size of each ensemble (N)
    except:
        print(
            "usage: python %s <database> <measurement> <exp_err and weight factors> <# of the ensemble>"
            % os.path.split(sys.argv[0])[1]
        )
        sys.exit(1)

    EXP_ERRO, WEIGHT_FACTOR = read_error_weights(error_wf)
    main(world, target)
    
    