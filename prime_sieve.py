from sympy import Integer, factorint, log, sieve, isprime
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import itertools
from rabinmiller import rabinmiller
# from aks import aks_test

# Generate composites from exponent vectors (memory-efficient representation)
def composite_from_exponents(primes, exponents):
    result = Integer(1)
    for p, e in zip(primes, exponents):
        result *= Integer(p)**e
    return result

# Generate nth composite explicitly with external indices
def generate_composites(primes, count, composites=None, indices=None):
    if composites is None:
        composites = [Integer(1)]
    if indices is None:
        indices = [0] * len(primes)

    while len(composites) < count:
        candidates = [primes[i] * composites[indices[i]] for i in range(len(primes))]
        next_comp = min(candidates)
        composites.append(next_comp)

        for i in range(len(primes)):
            if candidates[i] == next_comp:
                indices[i] += 1

    return composites, indices

# Improved logarithmic fit function (quadratic in log)
def log_fit(x, a, b, c):
    return a * np.log(x + 1)**2 + b * np.log(x + 1) + c

# Brute-force approach to find the nearest k-triple composite
def find_nearest_k_triple(primes, target):
    max_exponents = [int(np.log(int(target)) / np.log(p)) + 1 for p in primes]
    min_diff = float('inf')
    nearest = None

    for exponents in itertools.product(*[range(e + 1) for e in max_exponents]):
        composite = composite_from_exponents(primes, exponents)
        diff = abs(composite - target)
        if diff < min_diff:
            min_diff = diff
            nearest = composite

    return nearest

# Prime checking with Miller-Rabin and AKS for confirmation
def confirm_prime(candidate):
    if rabinmiller(int(candidate)):
        if isprime(candidate):
            return candidate
    return None

# Fo Sieve main function using initial known primes and outputting graph for each curve
def fo_sieve(prime_triple, initial_prime_count, num_primes_to_predict):
    known_primes = list(sieve.primerange(1, sieve[initial_prime_count]+1))
    composites, indices = generate_composites(prime_triple, known_primes[-1] + 2)

    indices_arr = np.arange(len(known_primes))
    logs_composites = np.log([float(composites[int(p)]) for p in known_primes])

    params, _ = curve_fit(log_fit, indices_arr, logs_composites)

    plt.figure()
    plt.plot(indices_arr, logs_composites, 'bo-', label='Known Composites Log Values')
    plt.plot(indices_arr, log_fit(indices_arr, *params), 'r--', label='Improved Logarithmic Fit')
    plt.title('Improved Curve Fit')
    plt.legend()
    plt.show()

    last_prime = known_primes[-1]
    predicted_primes = []
    current_index = len(known_primes)

    while len(predicted_primes) < num_primes_to_predict:
        predicted_log = log_fit(current_index, *params)
        predicted_composite = Integer(round(np.exp(predicted_log)))

        nearest_composite = find_nearest_k_triple(prime_triple, predicted_composite)

        while nearest_composite not in composites:
            composites, indices = generate_composites(prime_triple, len(composites) + 1, composites, indices)

        candidate = composites.index(nearest_composite) + 1
        if (candidate > last_prime and confirm_prime(candidate)):
            predicted_primes.append(candidate)
            last_prime = candidate
            current_index += 1
            continue

        confirmed = False
        done_left = False  
        print(f"candidate: {candidate}, last_prime: {last_prime}")
        delta = 1
        while not confirmed:
            left_candidate = candidate - delta
            right_candidate = candidate + delta if candidate + delta > last_prime else last_prime + delta
            done_left = (left_candidate <= last_prime)
            if (not done_left):
                print(f"left_candidate: {left_candidate}", end = ", ")
            print(f"right_candidate: {right_candidate}")
            if (not done_left and confirm_prime(left_candidate)):
                confirmed = True
                predicted_primes.append(left_candidate)
                last_prime = left_candidate
                print(f"Confirmed prime at index {current_index}: {left_candidate}")
            if (not confirmed and confirm_prime(right_candidate)):
                confirmed = True
                predicted_primes.append(right_candidate)
                last_prime = right_candidate
                print(f"Confirmed prime at index {current_index}: {right_candidate}")
            delta += 1

        current_index += 1

    return predicted_primes

# Example usage
prime_triple = [3, 5, 7]
initial_prime_count = 30
num_primes_to_predict = 20

predicted_primes = fo_sieve(prime_triple, initial_prime_count, num_primes_to_predict)
actual_primes = [x for x in sieve[31:51]]
print("Predicted primes:", predicted_primes)
print("Actual primes:", actual_primes)
