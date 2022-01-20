import hashlib
import string
import random
from Queue import Queue
from Pyro4 import expose


@expose
class Solver:
    def __init__(self, workers=None, input_file_name=None, output_file_name=None):
        self.input_file_name = input_file_name
        self.output_file_name = output_file_name
        self.workers = workers

    def solve(self):
        input_list = self.read_input()
        alphabet_tag = input_list.pop(0)
        alphabet = self.get_alphabet(alphabet_tag)
        pattern = input_list.pop(0)
        n = len(input_list)
        k = n // len(self.workers)
        j = n % len(self.workers)
        mapped = []
        for i in xrange(0, len(self.workers)):
            mapped.append(self.workers[i].boyer_moore(pattern, alphabet,
                                                      input_list[k * i:k * (i + 1) + j * ((i + 1) // len(self.workers))]))

        result = self.myreduce(mapped)
        self.write_output('\n'.join([pattern, str(result)]))

    @staticmethod
    def myreduce(mapped):
        output = 0
        for el in mapped:
            output += el.value
        return output

    @staticmethod
    @expose
    def boyer_moore(p, alphabet, texts):
        p_bm = BoyerMoore(p, alphabet)
        i = 0
        occurrences = 0
        for t in texts:
            while i < len(t) - len(p) + 1:
                shift = 1
                mismatched = False
                for j in range(len(p) - 1, -1, -1):
                    if p[j] != t[i + j]:
                        skip_bc = p_bm.bad_character_rule(j, t[i + j])
                        skip_gs = p_bm.good_suffix_rule(j)
                        shift = max(shift, skip_bc, skip_gs)
                        mismatched = True
                        break
                if not mismatched:
                    occurrences += 1
                    skip_gs = p_bm.match_skip()
                    shift = max(shift, skip_gs)
                i += shift
        return occurrences

    @staticmethod
    def get_alphabet(tags):
        tag_dict = {'d': string.digits,
                    'l': string.ascii_lowercase,
                    'u': string.ascii_uppercase,
                    'p': string.punctuation,
                    'w': string.whitespace}
        charset = ''.join([tag_dict[c] for c in tags if c in tag_dict])
        return charset if len(charset) > 0 else string.ascii_letters

    def read_input(self):
        input_list = []
        with open(self.input_file_name, 'r') as f:
            line = f.readline()
            while line:
                input_list.append(line.strip())
                line = f.readline()
        return input_list

    def write_output(self, output):
        with open(self.output_file_name, mode='w') as f:
            f.write(output)


def z_array(s):
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s) - 1)
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            z[1] += 1
        else:
            break
    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1
    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            for i in range(k, len(s)):
                if s[i] == s[i - k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                z[k] = zkp
            else:
                nmatch = 0
                for i in range(r + 1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    lp = [0] * len(p)
    for j in range(len(p) - 1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i - 1], lp[i])
    return l


def small_l_prime_array(n):
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i + 1:
            small_lp[len(n) - i - 1] = i + 1
    for i in range(len(n) - 2, -1, -1):
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i + 1]
    return small_lp


def good_suffix_table(p):
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i + 1
    return tab


class BoyerMoore(object):
    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        # Create map from alphabet characters to integers
        self.amap = {}
        for i in range(len(self.alphabet)):
            self.amap[self.alphabet[i]] = i
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        assert c in self.amap
        ci = self.amap[c]
        assert i > (self.bad_char[i][ci] - 1)
        return i - (self.bad_char[i][ci] - 1)

    def good_suffix_rule(self, i):
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        return len(self.small_l_prime) - self.small_l_prime[1]
