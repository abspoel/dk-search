import argparse
import re
from bisect import bisect_left
from collections import namedtuple
from math import ceil, log2

HighDkPrime = namedtuple('HighDkPrime', ['p', 'd0', 'd1', 'bits'])


def add_to_bests(bests, x, d0, d1):
    ix = 0
    if not bests:
        insert = True
    elif x < bests[0].p:
        insert = True
    else:
        ix = bisect_left([best.p for best in bests], x)
        if ix < len(bests) and bests[ix].p == x:
            return

        insert = True
        for best in bests[:ix]:
            if best.d0 >= d0 and best.d1 >= d1:
                insert = False

    if insert:
        bests[:] = [best for best in bests if best.p < x or best.d0 > d0 or best.d1 > d1]
        bests.insert(ix, HighDkPrime(x, d0, d1, int(ceil(log2(x)))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs='+')
    args = parser.parse_args()

    bests = list()
    for filename in args.files:
        lines = open(filename, "r").readlines()
        done = True
        k_values = None
        use_star = None

        for i, line in enumerate(lines):
            if line.startswith("=== Done"):
                break
        else:
            done = False
            regex = re.compile(r"(\d+) \((\d+) bits\), d(\d+)(\*?)=(\d+), d(\d+)(\*?)=(\d+)")
            for line in lines:
                if not line.strip():
                    continue
                match = regex.search(line)
                if not match:
                    continue

                f, k = match.group(3), match.group(6)
                star_match = match.group(4) == "*"
                if k_values is None:
                    k_values = f, k
                    use_star = star_match
                elif (f, k) != k_values or star_match != use_star:
                    raise ValueError("Line {} does not match k_values {} or use_star {} from before".format(line, k_values, use_star))

                p, d0, d1 = int(match.group(1)), int(match.group(5)), int(match.group(8))
                add_to_bests(bests, p, d0, d1)

        if done:
            regex = re.compile(r"^(\d+) \((\d+) bits\), d(\d+)(\*?)=(\d+), d(\d+)(\*?)=(\d+)")
            for line in lines[i+1:]:
                if not line.strip():
                    continue
                match = regex.match(line)
                if not match:
                    raise ValueError("Can't parse line {}".format(line))
                f, k = match.group(3), match.group(6)
                star_match = match.group(4) == "*"
                if k_values is None:
                    k_values = f, k
                    use_star = star_match
                elif (f, k) != k_values or star_match != use_star:
                    raise ValueError("Line {} does not match k_values {} or use_star {} from before".format(line, k_values, use_star))

                p, d0, d1 = int(match.group(1)), int(match.group(5)), int(match.group(8))
                add_to_bests(bests, p, d0, d1)

    star = "*" if use_star else ""
    d1_best = None
    d1_bests = []
    for best in bests:
        print("{} ({} bits), d{}{}={}, d{}{}={}".format(best.p, best.bits, k_values[0], star, best.d0, k_values[1], star, best.d1))
        if d1_best is None:
            d1_best = best
        elif best.d1 > d1_best.d1:
            #if best.bits > d1_best.bits:
            d1_bests.append(d1_best)
            d1_best = best
    d1_bests.append(d1_best)
    
    print("\nBest d{}:".format(k_values[1]))
    for best in d1_bests:
        print("{} ({} bits), d{}{}={}, d{}{}={}".format(best.p, best.bits, k_values[0], star, best.d0, k_values[1], star, best.d1))



if __name__ == "__main__":
    main()
