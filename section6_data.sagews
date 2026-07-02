from sage.all import *
from math import ceil
from collections import OrderedDict
import re
import sys

def has_rational_root_of_division_polynomial(E, k):
    R.<x> = PolynomialRing(QQ)
    psik = R(E.division_polynomial(k))
    return any(factor.degree() == 1 for factor, exponent in psik.factor())


def P_value(E, k, m, skip_bad_primes=True):
    NE = E.conductor()
    count = 0

    for p in prime_range(m + 1):
        if skip_bad_primes and NE % p == 0:
            continue

        ap = E.ap(p)
        Ap = (p + 1)^2 - ap^2

        if gcd(Ap, k) == 1:
            count += 1

    return count


def print_in_columns(items, ncols=7, width=24):
    if len(items) == 0:
        print("No curves found.")
        return

    nrows = ceil(len(items) / ncols)

    for i in range(nrows):
        row = ""
        for j in range(ncols):
            idx = i + j*nrows
            if idx < len(items):
                row += "{:<{w}}".format(items[idx], w=width)
        print(row)


def isogeny_class(label):
    return re.match(r"(\d+[a-z]+)", label).group(1)


def print_grouped_results(results, ncols=7, width=24):

    if len(results) == 0:
        print("No curves found.")
        return

    groups = OrderedDict()

    for label, count, answer in results:
        key = isogeny_class(label)

        if key not in groups:
            groups[key] = []

        groups[key].append(f"{label} ({count}, {answer})")

    first = True

    for items in groups.values():
        if not first:
            print("-" * 80)

        print_in_columns(items, ncols=ncols, width=width)
        first = False


def search_curves(N, k, m, threshold=10,
                  skip_bad_primes=True, ncols=7, width=24):
    """
    Yes = k-division polynomial has a rational root.
    No  = k-division polynomial has no rational root.
    """
    results = []

    for E in cremona_curves(range(1, N + 1)):
        label = E.cremona_label()
        count = P_value(E, k, m, skip_bad_primes=skip_bad_primes)

        if count < threshold:
            has_root = has_rational_root_of_division_polynomial(E, k)
            answer = "yes" if has_root else "NO"

            results.append((label, count, answer))

    print("Entry format: Cremona label (Permutation primes, Rational k-division root?)")
    print("Permutation primes = # of primes p <=", m,
          "with gcd((p+1)^2-a_p^2,", k, ") = 1.")
    print("Yes = k-th division polynomial has a rational root, No = none")
    print()

    print_grouped_results(results, ncols=ncols, width=width)

    return results


def run_for_k_range(r, N=1000, m=1000, threshold=10,
                    skip_bad_primes=True, ncols=7, width=24):
    all_results = {}

    for k in range(2, r + 1):
        print("\n" + "=" * 80)
        print("k =", k)
        print("=" * 80)

        results_k = search_curves(
            N=N,
            k=k,
            m=m,
            threshold=threshold,
            skip_bad_primes=skip_bad_primes,
            ncols=ncols,
            width=width
        )

        all_results[k] = results_k

    return all_results
