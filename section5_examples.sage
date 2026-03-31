from sage.all import *
from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective

# Usage:
#   sage section5_examples.sage
# In a restricted sandbox, use:
#   DOT_SAGE=/tmp/.sage sage section5_examples.sage


OUTPUT_WIDTH = 78


def banner(title):
    line = "=" * OUTPUT_WIDTH
    print(line)
    print(title)
    print(line)


def subbanner(title):
    line = "-" * OUTPUT_WIDTH
    print("\n" + line)
    print(title)
    print(line)


def poly_ring(base_ring):
    return PolynomialRing(base_ring, "x")


def lattes_map_polynomials(E, k, base_ring=QQ):
    Rx = poly_ring(base_ring)
    x = Rx.gen()
    E0 = E.change_ring(base_ring)
    num = E0._multiple_x_numerator(k, x=x)
    den = E0._multiple_x_denominator(k, x=x)
    g = num.gcd(den)
    if g.degree() > 0:
        num = num // g
        den = den // g
    return num, den


def rational_function_string(num, den):
    return f"({num}) / ({den})"


def homogenize(poly, X, Y, degree):
    out = X.parent().zero()
    for exponent, coeff in poly.dict().items():
        out += coeff * X**exponent * Y**(degree - exponent)
    return out


def lattes_system_mod_p(E, k, p):
    F = GF(p)
    P1 = ProjectiveSpace(F, 1, names=("X", "Y"))
    X, Y = P1.coordinate_ring().gens()
    num, den = lattes_map_polynomials(E, k, F)
    degree = max(num.degree(), den.degree())
    num_h = homogenize(num, X, Y, degree)
    den_h = homogenize(den, X, Y, degree)
    return DynamicalSystem_projective([num_h, den_h], domain=P1)


def point_label(point):
    if point[1] == 0:
        return "infinity"
    return ZZ(point[0] / point[1])


def p1_value_table(E, k, p):
    F = GF(p)
    phi = lattes_system_mod_p(E, k, p)
    P1 = phi.domain()
    xs = list(range(p)) + ["infinity"]
    values = [point_label(phi(P1([F(i), 1]))) for i in range(p)]
    values.append(point_label(phi(P1([1, 0]))))
    return xs, values


def permutation_data(E, D, k, primes):
    rows = []
    for p in primes:
        ap = ZZ(E.ap(p))
        gcd_value = ZZ(gcd((p + 1) ** 2 - ap**2, k))
        _, values = p1_value_table(E, k, p)
        image_size = len(set(values))
        permutes = image_size == p + 1
        criterion = gcd_value == 1
        if permutes != criterion:
            raise RuntimeError(
                f"Criterion mismatch for D={D}, k={k}, p={p}: "
                f"gcd={gcd_value}, image_size={image_size}"
            )
        rows.append(
            (
                p,
                ZZ(kronecker(D, p)),
                ap,
                gcd_value,
                image_size,
                "Yes" if permutes else "No",
            )
        )
    return rows


def print_rows(headers, rows):
    widths = [len(str(h)) for h in headers]
    for row in rows:
        for i, entry in enumerate(row):
            widths[i] = max(widths[i], len(str(entry)))
    print("  ".join(str(h).ljust(widths[i]) for i, h in enumerate(headers)))
    print("  ".join("-" * widths[i] for i in range(len(headers))))
    for row in rows:
        print("  ".join(str(entry).ljust(widths[i]) for i, entry in enumerate(row)))


def torsion_label(E):
    invariants = tuple(ZZ(n) for n in E.torsion_subgroup().invariants())
    if not invariants:
        return "C1"
    return " x ".join(f"C{n}" for n in invariants)


def torsion_points_with_orders(E):
    rows = []
    for P in E.torsion_points():
        rows.append((str(P), ZZ(P.order())))
    rows.sort(key=lambda item: (item[1], item[0]))
    return rows


def rational_roots(poly):
    return sorted(poly.roots(QQ, multiplicities=False))


def rational_lifts_at_x(E, x0):
    points = []
    for P in E.lift_x(QQ(x0), all=True):
        points.append((str(P), ZZ(P.order())))
    points.sort(key=lambda item: (item[1], item[0]))
    return points


def quadratic_factor_discriminants(poly):
    data = []
    for factor, multiplicity in poly.factor():
        if factor.degree() == 2:
            data.append((factor, multiplicity, factor.discriminant()))
    data.sort(key=lambda item: str(item[0]))
    return data


def fraction_field_function(num, den):
    Rx = poly_ring(QQ)
    K = Frac(Rx)
    return K(num) / K(den)


def report_summary(cases):
    subbanner("Summary of Curves Used in Section 5")
    rows = []
    for case in cases:
        E = case["curve"]
        rows.append(
            (
                case["D"],
                case["curve_short"],
                E.cm_discriminant(),
                E.j_invariant(),
                torsion_label(E),
            )
        )
    print_rows(["D", "curve", "cm_disc", "j", "torsion"], rows)


def report_division_polynomial(E, n):
    subbanner(f"Division polynomial psi_{n}(x)")
    poly = E.division_polynomial(n)
    print(f"psi_{n}(x) = {poly}")
    print(f"factorization over QQ = {poly.factor()}")
    print(f"irreducible over QQ? {poly.is_irreducible()}")
    roots = rational_roots(poly)
    print(f"rational roots = {roots}")
    if roots:
        print("rational lifts on the chosen Weierstrass model:")
        for root in roots:
            lifts = rational_lifts_at_x(E, root)
            if lifts:
                print(f"  x = {root}: {lifts}")
            else:
                print(
                    f"  x = {root}: no rational lift on E(Q); "
                    f"this is still a rational x-coordinate of n-torsion."
                )
    discs = quadratic_factor_discriminants(poly)
    if discs:
        print("quadratic factors and discriminants:")
        for factor, multiplicity, disc in discs:
            print(f"  factor = {factor}, multiplicity = {multiplicity}, discriminant = {disc}")


def report_lattes_map(E, k, print_formula):
    subbanner(f"L_{k} for the x-coordinate map")
    num, den = lattes_map_polynomials(E, k)
    print(f"numerator degree = {num.degree()}")
    print(f"denominator degree = {den.degree()}")
    if print_formula:
        print(f"L_{k}(x) = {rational_function_string(num, den)}")


def report_prime_table(E, D, k, primes, title):
    subbanner(title)
    rows = permutation_data(E, D, k, primes)
    print_rows(["p", "symbol", "a_p", "gcd", "image_size", "permutes?"], rows)
    print("criterion check: in every displayed row, gcd((p+1)^2-a_p^2, k)=1")
    print("if and only if the direct action on P^1(F_p) is a permutation.")


def report_value_table(E, k, p, title):
    subbanner(title)
    xs, values = p1_value_table(E, k, p)
    print_rows(["entry"] + [str(x) for x in xs], [("x",) + tuple(xs), (f"L_{k}(x)",) + tuple(values)])
    print(f"image size = {len(set(values))} out of {p + 1}")


def report_composition(E, k1, k2, k12):
    subbanner(f"Composition check for L_{k12}")
    num1, den1 = lattes_map_polynomials(E, k1)
    num2, den2 = lattes_map_polynomials(E, k2)
    num12, den12 = lattes_map_polynomials(E, k12)
    f1 = fraction_field_function(num1, den1)
    f2 = fraction_field_function(num2, den2)
    f12 = fraction_field_function(num12, den12)
    print(f"L_{k1}(L_{k2}(x)) == L_{k12}(x)? {f1(f2) == f12}")
    print(f"L_{k2}(L_{k1}(x)) == L_{k12}(x)? {f2(f1) == f12}")


def report_torsion_points(E):
    subbanner("Rational torsion points")
    rows = torsion_points_with_orders(E)
    print_rows(["point", "order"], rows)


cases = [
    {
        "title": "Case D in {-4, -16}",
        "D": -4,
        "curve_short": "y^2=x^3+x",
        "curve": EllipticCurve([0, 0, 0, 1, 0]),
        "division_jobs": [3],
        "lattes_jobs": [(3, True)],
        "prime_tables": [("Permutation data for L_3 at selected primes", 3, [5, 7, 11, 13, 19, 23, 29, 31])],
        "value_tables": [("Values of L_3 on P^1(F_7)", 3, 7)],
    },
    {
        "title": "Case D = -8",
        "D": -8,
        "curve_short": "y^2=x^3+x^2-3x+1",
        "curve": EllipticCurve([0, 1, 0, -3, 1]),
        "division_jobs": [],
        "lattes_jobs": [(3, True)],
        "prime_tables": [],
        "value_tables": [("Values of L_3 on P^1(F_13)", 3, 13)],
    },
    {
        "title": "Case D in {-7, -28}",
        "D": -7,
        "curve_short": "y^2=x^3-35x+98",
        "curve": EllipticCurve([0, 0, 0, -35, 98]),
        "division_jobs": [3],
        "lattes_jobs": [(3, True)],
        "prime_tables": [("Permutation data for L_3 at selected primes", 3, [5, 11, 13, 17, 19, 23, 29, 31])],
        "value_tables": [],
    },
    {
        "title": "Case D in {-19, -43, -67, -163}",
        "D": -19,
        "curve_short": "y^2+y=x^3-38x+90",
        "curve": EllipticCurve([0, 0, 1, -38, 90]),
        "division_jobs": [],
        "lattes_jobs": [(3, True)],
        "prime_tables": [("Permutation data for L_3 at selected primes", 3, [5, 7, 11, 13, 17, 23, 29, 41])],
        "value_tables": [("Values of L_3 on P^1(F_5)", 3, 5)],
    },
    {
        "title": "Case D = -12",
        "D": -12,
        "curve_short": "y^2=x^3-15x+22",
        "curve": EllipticCurve([0, 0, 0, -15, 22]),
        "division_jobs": [],
        "lattes_jobs": [(5, True)],
        "prime_tables": [],
        "value_tables": [("Values of L_5 on P^1(F_11)", 5, 11)],
        "show_torsion_points": True,
    },
    {
        "title": "Case D = -27",
        "D": -27,
        "curve_short": "y^2=x^3-120x+506",
        "curve": EllipticCurve([0, 0, 0, -120, 506]),
        "division_jobs": [2, 3],
        "lattes_jobs": [(5, True), (2, True)],
        "prime_tables": [
            ("Permutation data for L_5 at selected primes", 5, [5, 11, 17, 23, 29, 41, 47]),
            ("Permutation data for L_2 at selected primes", 2, [7, 13, 19, 31, 37, 43, 61]),
        ],
        "value_tables": [],
    },
    {
        "title": "Case D = -3",
        "D": -3,
        "curve_short": "y^2=x^3+2",
        "curve": EllipticCurve([0, 0, 0, 0, 2]),
        "division_jobs": [2, 3],
        "lattes_jobs": [(2, True)],
        "prime_tables": [("Permutation data for L_2 at selected primes", 2, [5, 7, 11, 13, 19, 31, 37, 43])],
        "value_tables": [("Values of L_2 on P^1(F_7)", 2, 7)],
    },
    {
        "title": "Counterexample case D = -11",
        "D": -11,
        "curve_short": "y^2=x^3-264x+1694",
        "curve": EllipticCurve([0, 0, 0, -264, 1694]),
        "division_jobs": [2, 3],
        "lattes_jobs": [(2, True), (3, True), (6, True), (7, True)],
        "prime_tables": [
            ("Permutation data for L_6 at selected primes", 6, [5, 7, 13, 17, 23, 31, 41, 47]),
            ("Permutation data for L_7 at selected primes", 7, [5, 7, 13, 17, 23, 31]),
        ],
        "value_tables": [
            ("Values of L_6 on P^1(F_5)", 6, 5),
            ("Values of L_6 on P^1(F_7)", 6, 7),
            ("Values of L_7 on P^1(F_5)", 7, 5),
        ],
        "show_composition": (2, 3, 6),
    },
]


banner("Section 5 computation file for main_11_03_2026/main.tex")
print("All formulas and tables below are computed directly in SageMath from the")
print("elliptic curves used in Section 5. No manuscript table entries are fed")
print("back into the computation.")

report_summary(cases)

for case in cases:
    E = case["curve"]
    subbanner(case["title"])
    print(f"curve = {case['curve_short']}")
    print(f"Sage model = {E}")
    print(f"CM discriminant = {E.cm_discriminant()}")
    print(f"j-invariant = {E.j_invariant()}")
    print(f"torsion subgroup = {torsion_label(E)}")

    if case.get("show_torsion_points"):
        report_torsion_points(E)

    for n in case["division_jobs"]:
        report_division_polynomial(E, n)

    for k, print_formula in case["lattes_jobs"]:
        report_lattes_map(E, k, print_formula)

    if "show_composition" in case:
        report_composition(E, *case["show_composition"])

    for title, k, primes in case["prime_tables"]:
        report_prime_table(E, case["D"], k, primes, title)

    for title, k, p in case["value_tables"]:
        report_value_table(E, k, p, title)


banner("End of Section 5 computation")
