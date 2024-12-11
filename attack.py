import os
import sys
import time
import logging

from sage.all import *
from sage.crypto.util import *

DEBUG_ROOTS = None
BOUND_CHECK = False
USE_FLATTER = True
ACLOG_CLEAR = True

log_file = 'attack.log'  
if ACLOG_CLEAR and os.path.exists(log_file):  
    os.remove(log_file) 
logger = logging.getLogger(__name__)
logging.basicConfig(filename = log_file, level = logging.DEBUG, format = '%(asctime)s - %(levelname)s - %(message)s')


def create_lattice(pr, shifts, bounds, order="invlex", sort_shifts_reverse=False, sort_monomials_reverse=False):
    """
    Creates a lattice from a list of shift polynomials.
    :param pr: the polynomial ring
    :param shifts: the shifts
    :param bounds: the bounds
    :param order: the order to sort the shifts/monomials by
    :param sort_shifts_reverse: set to true to sort the shifts in reverse order
    :param sort_monomials_reverse: set to true to sort the monomials in reverse order
    :return: a tuple of lattice and list of monomials
    """
    logging.debug(f"Creating a lattice with {len(shifts)} shifts ({order = }, {sort_shifts_reverse = }, {sort_monomials_reverse = })...")
    if pr.ngens() > 1:
        pr_ = pr.change_ring(ZZ, order=order)
        shifts = [pr_(shift) for shift in shifts]

    monomials = set()
    for shift in shifts:
        monomials.update(shift.monomials())

    shifts.sort(reverse=sort_shifts_reverse)
    monomials = sorted(monomials, reverse=sort_monomials_reverse)
    L = matrix(ZZ, len(shifts), len(monomials))
    for row, shift in enumerate(shifts):
        for col, monomial in enumerate(monomials):
            L[row, col] = shift.monomial_coefficient(monomial) * monomial(*bounds)

    monomials = [pr(monomial) for monomial in monomials]
    return L, monomials


def reduce_lattice(L, delta=0.8):
    """
    Reduces a lattice basis using a lattice reduction algorithm (currently LLL).
    :param L: the lattice basis
    :param delta: the delta parameter for LLL (default: 0.8)
    :return: the reduced basis
    """
    # logging.debug(f"Reducing a {L.nrows()} x {L.ncols()} lattice...")
    # return L.LLL(delta)
    start_time = time.perf_counter()
    if USE_FLATTER:
        from subprocess import check_output
        from re import findall
        LL = "[[" + "]\n[".join(" ".join(map(str, row)) for row in L) + "]]"
        ret = check_output(["flatter"], input = LL.encode())
        L_reduced = matrix(L.nrows(), L.ncols(), map(int, findall(rb"-?\d+", ret)))
    else:
        L_reduced = L.LLL(delta)
    end_time = time.perf_counter()
    reduced_time = end_time - start_time
    logging.info(f"Reducing a {L.nrows()} x {L.ncols()} lattice within {reduced_time:.3f} seconds...")
    return L_reduced


def reconstruct_polynomials(B, f, modulus, monomials, bounds, preprocess_polynomial=lambda x: x, divide_gcd=True):
    """
    Reconstructs polynomials from the lattice basis in the monomials.
    :param B: the lattice basis
    :param f: the original polynomial (if set to None, polynomials will not be divided by f if possible)
    :param modulus: the original modulus
    :param monomials: the monomials
    :param bounds: the bounds
    :param preprocess_polynomial: a function which preprocesses a polynomial before it is added to the list (default: identity function)
    :param divide_gcd: if set to True, polynomials will be pairwise divided by their gcd if possible (default: True)
    :return: a list of polynomials
    """
    divide_original = f is not None
    modulus_bound = modulus is not None
    logging.debug(f"Reconstructing polynomials ({divide_original = }, {modulus_bound = }, {divide_gcd = })...")
    polynomials = []
    for row in range(B.nrows()):
        norm_squared = 0
        w = 0
        polynomial = 0
        for col, monomial in enumerate(monomials):
            if B[row, col] == 0:
                continue
            norm_squared += B[row, col] ** 2
            w += 1
            assert B[row, col] % monomial(*bounds) == 0
            polynomial += B[row, col] * monomial // monomial(*bounds)

        # Equivalent to norm >= modulus / sqrt(w)
        # Use BOUND_CHECK = False to achieve a successful attack
        if BOUND_CHECK and modulus_bound and norm_squared * w >= modulus ** 2:
            logging.debug(f"Row {row} is too large, ignoring...")
            continue

        polynomial = preprocess_polynomial(polynomial)

        if divide_original and polynomial % f == 0:
            logging.debug(f"Original polynomial divides reconstructed polynomial at row {row}, dividing...")
            polynomial //= f

        if divide_gcd:
            for i in range(len(polynomials)):
                g = gcd(polynomial, polynomials[i])
                # TODO: why are we only allowed to divide out g if it is constant?
                if g != 1 and g.is_constant():
                    logging.debug(f"Reconstructed polynomial has gcd {g} with polynomial at {i}, dividing...")
                    polynomial //= g
                    polynomials[i] //= g

        if polynomial.is_constant():
            logging.debug(f"Polynomial at row {row} is constant, ignoring...")
            continue

        if DEBUG_ROOTS is not None:
            logging.debug(f"Polynomial at row {row} roots check: {polynomial(*DEBUG_ROOTS)}")

        polynomials.append(polynomial)

    logging.debug(f"Reconstructed {len(polynomials)} polynomials")
    return polynomials


def find_roots_univariate(x, polynomial):
    """
    Returns a generator generating all roots of a univariate polynomial in an unknown.
    :param x: the unknown
    :param polynomial: the polynomial
    :return: a generator generating dicts of (x: root) entries
    """
    if polynomial.is_constant():
        return

    for root in polynomial.roots(multiplicities=False):
        if root != 0:
            yield {x: int(root)}


def find_roots_gcd(pr, polynomials):
    """
    Returns a generator generating all roots of a polynomial in some unknowns.
    Uses pairwise gcds to find trivial roots.
    :param pr: the polynomial ring
    :param polynomials: the reconstructed polynomials
    :return: a generator generating dicts of (x0: x0root, x1: x1root, ...) entries
    """
    if pr.ngens() != 2:
        return

    logging.debug("Computing pairwise gcds to find trivial roots...")
    x, y = pr.gens()
    for i in range(len(polynomials)):
        for j in range(i):
            g = gcd(polynomials[i], polynomials[j])
            if g.degree() == 1 and g.nvariables() == 2 and g.constant_coefficient() == 0:
                # g = ax + by
                a = int(g.monomial_coefficient(x))
                b = int(g.monomial_coefficient(y))
                yield {x: b, y: a}
                yield {x: -b, y: a}


def find_roots_groebner(pr, polynomials):
    """
    Returns a generator generating all roots of a polynomial in some unknowns.
    Uses Groebner bases to find the roots.
    :param pr: the polynomial ring
    :param polynomials: the reconstructed polynomials
    :return: a generator generating dicts of (x0: x0root, x1: x1root, ...) entries
    """
    # We need to change the ring to QQ because groebner_basis is much faster over a field.
    # We also need to change the term order to lexicographic to allow for elimination.
    gens = pr.gens()
    s = Sequence(polynomials, pr.change_ring(QQ, order="lex"))
    while len(s) > 0:
        G = s.groebner_basis()
        logging.debug(f"Sequence length: {len(s)}, Groebner basis length: {len(G)}")
        if len(G) == len(gens):
            logging.debug(f"Found Groebner basis with length {len(gens)}, trying to find roots...")
            roots = {}
            for polynomial in G:
                vars = polynomial.variables()
                if len(vars) == 1:
                    for root in find_roots_univariate(vars[0], polynomial.univariate_polynomial()):
                        roots |= root

            if len(roots) == pr.ngens():
                yield roots
                return

            logging.debug(f"System is underdetermined, trying to find constant root...")
            G = Sequence(s, pr.change_ring(ZZ, order="lex")).groebner_basis()
            vars = tuple(map(lambda x: var(x), gens))
            for solution_dict in solve([polynomial(*vars) for polynomial in G], vars, solution_dict=True):
                logging.debug(solution_dict)
                found = False
                roots = {}
                for i, v in enumerate(vars):
                    s = solution_dict[v]
                    if s.is_constant():
                        if not s.is_zero():
                            found = True
                        roots[gens[i]] = int(s) if s.is_integer() else int(s) + 1
                    else:
                        roots[gens[i]] = 0
                if found:
                    yield roots
                    return

            return
        else:
            # Remove last element (the biggest vector) and try again.
            s.pop()


def find_roots_resultants(gens, polynomials):
    """
    Returns a generator generating all roots of a polynomial in some unknowns.
    Recursively computes resultants to find the roots.
    :param polynomials: the reconstructed polynomials
    :param gens: the unknowns
    :return: a generator generating dicts of (x0: x0root, x1: x1root, ...) entries
    """
    if len(polynomials) == 0:
        return

    if len(gens) == 1:
        if polynomials[0].is_univariate():
            yield from find_roots_univariate(gens[0], polynomials[0].univariate_polynomial())
    else:
        resultants = [polynomials[0].resultant(polynomials[i], gens[0]) for i in range(1, len(gens))]
        for roots in find_roots_resultants(gens[1:], resultants):
            for polynomial in polynomials:
                polynomial = polynomial.subs(roots)
                if polynomial.is_univariate():
                    for root in find_roots_univariate(gens[0], polynomial.univariate_polynomial()):
                        # Show a root 
                        logging.debug(f"Now root is {root}")
                        yield roots | root


def find_roots_variety(pr, polynomials):
    """
    Returns a generator generating all roots of a polynomial in some unknowns.
    Uses the Sage variety (triangular decomposition) method to find the roots.
    :param pr: the polynomial ring
    :param polynomials: the reconstructed polynomials
    :return: a generator generating dicts of (x0: x0root, x1: x1root, ...) entries
    """
    # We need to change the ring to QQ because variety requires a field.
    s = Sequence([], pr.change_ring(QQ))
    # We use more polynomials (i.e., poly_number) to find the roots, we can further tweak it
    poly_number = int(len(polynomials) * 0.5)
    for i in range(poly_number):
        s.append(polynomials[i])
    I = s.ideal()
    dim = I.dimension()
    logging.debug(f"{I.groebner_basis() = }")
    logging.debug(f"Sequence length: {len(s)}, Ideal dimension: {dim}")
    if dim == 0:
        logging.debug("Found ideal with dimension 0, computing variety...")
        logging.debug(f"The variety is {I.variety(ring=ZZ)}")
        for roots in I.variety(ring=ZZ):
            yield {k: int(v) for k, v in roots.items()}
        return
    elif dim == 1:
        logging.debug("Found ideal with dimension 1...")
        logging.debug(f"{I.groebner_basis()}")
        yield I.groebner_basis()[0]
        return 


def find_roots(pr, polynomials, method="groebner"):
    """
    Returns a generator generating all roots of a polynomial in some unknowns.
    The method used depends on the method parameter.
    :param pr: the polynomial ring
    :param polynomials: the reconstructed polynomials
    :param method: the method to use, can be "groebner", "resultants", or "variety" (default: "groebner")
    :return: a generator generating dicts of (x0: x0root, x1: x1root, ...) entries
    """
    if pr.ngens() == 1:
        logging.debug("Using univariate polynomial to find roots...")
        for polynomial in polynomials:
            yield from find_roots_univariate(pr.gen(), polynomial)
    else:
        # Always try this method because it can find roots the others can't.
        yield from find_roots_gcd(pr, polynomials)

        if method == "groebner":
            logging.debug("Using Groebner basis method to find roots...")
            yield from find_roots_groebner(pr, polynomials)
        elif method == "resultants":
            logging.debug("Using resultants method to find roots...")
            yield from find_roots_resultants(pr.gens(), polynomials)
        elif method == "variety":
            logging.debug("Using variety method to find roots...")
            yield from find_roots_variety(pr, polynomials)


def modular_bivariate(N, e, delta, s, desired_solution, roots_method="resultants"):
    """
    Computes small modular roots of a bivariate polynomial.
    More information: Zheng M., Kunihiro N., Yao Y., "Cryptanalysis of the RSA variant based on cubic Pell equation"
    :param N: the modulus
    :param e: the public exponent
    :param delta: the given size on the private exponent
    :param s: the parameter s
    :param desired_solution: used for debug info
    :param roots_method: the suitable method to use to find roots (default: "resultants")
    :return: a generator generating small roots (tuples of x, y roots) of the polynomial
    """
    fr = ZZ["x", "y"]
    x, y = fr.gens()
    f = x * y ** 2 + (N + 1) * x * y + (N ** 2 - N + 1) * x + 1
    f = f.change_ring(ZZ)
    logging.info(f"Polynomial to be solved: {f = }")
    pr = ZZ["x", "y", "z"]
    x, y, z = pr.gens()
    qr = pr.quotient(1 + x * y ** 2 - z)
    X = int(e * N ** delta / (N ** 2))
    Y = int(3 * N ** (1 / 2))
    Z = X * Y ** 2 + 1

    tau = min(2 - 2 * delta, 1)
    logging.info(f"Trying {s = }, {tau = :.3f}...")
    logging.info("Generating shifts...")
    shifts = []
    for i in range(s + 1):
        for k in range(s - i + 1):
            for j in [0, 1]:
                logging.debug(f"indices: {i = }, {j = }, {k= }")
                g = x ** i * y ** j * f ** k * e ** (s - k)  
                g = qr(g).lift()
                shifts.append(g)
    i = 0
    for k in range(s + 1):
        for j in range(2, floor(tau * k) + 2):
            logging.debug(f"indices: {i = }, {j = }, {k= }")
            h = x ** i * y ** j * f ** k * e ** (s - k)  
            h = qr(h).lift()
            shifts.append(h)
    
    monomials = set()
    for shift in shifts:
        monomials.add(shift.lm())
    logging.debug(f"The monomials: {monomials}")
    x0, y0 = desired_solution
    z0 = x0 * y0 ** 2 + 1
    for idx, shift in enumerate(shifts):
        logging.debug(f"Test for {idx + 1}: {shift(x0, y0, z0) % (e ** s)= }")

    logging.info("Generating the lattice...")
    L, monomials = create_lattice(pr, shifts, [X, Y, Z])
    logging.info("Reducing the lattice...")
    L = reduce_lattice(L)
    logging.debug(f"Test for original {f(x0, y0) % e = }")
    fr = f.parent()
    x, y = fr.gens()
    polynomials = reconstruct_polynomials(L, f, e ** s, monomials, [X, Y, Z], preprocess_polynomial = lambda p: p(x, y, 1 + x * y ** 2))
    for idx, poly in enumerate(polynomials):
        result = "0" if poly(x0, y0) % (e ** s) == 0 else "!@#$%"
        logging.debug(f"Test for {idx + 1} reconstructed poly(x0, y0) % (e ** s) = {result}")
        result = "0" if poly(x0, y0) == 0 else "!@#$%"
        logging.debug(f"Test for {idx + 1} reconstructed poly(x0, y0) = {result} over the integers")
    start_time = time.perf_counter()
    solutions = list(find_roots(fr, polynomials, method=roots_method))
    end_time = time.perf_counter()
    solution_time = end_time - start_time
    logging.debug(f"Finding roots within {solution_time:.3f} seconds...")
    for xy in solutions:
        x0 = xy[x]
        y0 = xy[y]
        if x0 != 0 and y0 != 0:
            logging.info(f"Found one possible solution: {x0 = } and {y0 = }")
            return x0, y0

    return None


def generate_cPRSA_instance(modulus_bit_length, delta, max_attempts = 10):
    """
    Generate an cubic Pell RSA instance with given bit-lengths of the modulus and private key d
    :param modulus_bit_length: the bit length of the modulus
    :param delta: a given size on the private exponent (d is roughly N^delta)
    :param max_attempts: the maximum number of attempts to generate cubic Pell RSA instance. (default: 10)
    :return: a tuple containing the public key (N, e)
    """
    e = d = N = Integer(1)
    attempts = 0
    d_bit_length = ceil(modulus_bit_length * delta)
    logging.info(f"Generating RSA instance with {modulus_bit_length}-bit modulus and {d_bit_length}-bit private key d...")
    while attempts < max_attempts and N.nbits() != modulus_bit_length and d.nbits() != d_bit_length:
        attempts += 1
        set_random_seed(int(time.time()))
        prime_bit_length = modulus_bit_length // 2
        p = random_blum_prime(2 ** (prime_bit_length - 1), 2 ** prime_bit_length - 1)
        q = random_blum_prime(2 ** (prime_bit_length - 1), 2 ** prime_bit_length - 1)
        N = p * q
        phi = (p ** 2 + p + 1) * (q ** 2 + q + 1)
        d = random_blum_prime(2 ** (d_bit_length - 1), 2 ** d_bit_length - 1)
        e = inverse_mod(d, phi)
        k = (e * d - 1) / phi
        p_q = p + q
        cPRSA_instance = [N, e, d, k, p, q] 
        desired_solution = [k, p_q]
        logging.debug(f'N: {N}')
        logging.debug(f'e: {e}')
        logging.debug(f'd: {d}')
        logging.debug(f'k: {k}')
        logging.debug(f'p: {p}')
        logging.debug(f'q: {q}')
        return cPRSA_instance, desired_solution
    logging.warning(f"Failed to generate cubic Pell RSA instance after {max_attempts} attempts...")
    return None


def attack_cPRSA_instance(modulus_bit_length, delta, s = 3):
    """
    Small private key attack on cubic Pell RSA instance with given parameters
    :param modulus_bit_length: the bit length of the modulus.
    :param delta: the ratio of the bit length of the private key to the modulus bit length.
    :param s: the given parameter for controlling the lattice dimension. (default: 3)
    :return: 1 if attack succeeds else 0
    """
    result = generate_cPRSA_instance(modulus_bit_length, delta)
    if result is not None:
        cPRSA_instance, desired_solution = result
    else:
        print(f"Sorry, cannot generate a cubic Pell RSA instance...")
        return 0, 0
    N, e = cPRSA_instance[0], cPRSA_instance[1]
    print(f"The parameters:\n{N = }\n{e = }")

    start_time = time.perf_counter()
    solution = modular_bivariate(N, e, delta, s, desired_solution)
    end_time = time.perf_counter()
    test_time = end_time - start_time
    if solution is not None:
        _, y0 = solution
        z = sqrt(y0 ** 2 - 4 * N)
        p = (y0 + z) // 2
        q = (y0 - z) // 2
        if p * q == N:
            logging.info(f"Succeeded!")
            logging.info(f"Found p = {p}")
            logging.info(f"Found q = {q}")
            print(f"Found primes:\n{p}\n{q}")
            return 1, test_time
        else:
            logging.info(f"Failed!")
            return 0, test_time
    else:
        print(f"Sorry, cannot attack this cubic Pell RSA instance...")
        return 0, test_time


if __name__ == "__main__":

    assert len(sys.argv) == 4, f"Usage: sage -python attack.py <modulus_bit_length> <delta> <s>"

    modulus_bit_length, delta, s = int(sys.argv[1]), RR(sys.argv[2]), int(sys.argv[3])
    result, test_time = attack_cPRSA_instance(modulus_bit_length, delta, s)
    if result:
        print(f"The attack costs {test_time:.3f} seconds...")