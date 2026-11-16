import math
import random
import statistics

K = 0.25
T = list(range(7))
C = [round(math.exp(-K * t), 3) for t in T]
I = 3  
EXACT = -K * math.exp(-K * T[I])


def forward(c, i, h=1):
    return (c[i + h] - c[i]) / h


def backward(c, i, h=1):
    return (c[i] - c[i - h]) / h


def central(c, i, h=1):
    return (c[i + h] - c[i - h]) / (2 * h)


def report(h):
    print(f"\nStep h = {h} min   (exact dC/dt = {EXACT:.4f} mol/L/min)")
    print(f"{'method':10s}{'estimate':>10s}{'error':>10s}{'% error':>10s}{'k = -rate/C':>14s}")
    errors = {}
    for name, fn in [("forward", forward), ("backward", backward), ("central", central)]:
        d = fn(C, I, h)
        err = d - EXACT
        errors[name] = err
        print(f"{name:10s}{d:10.4f}{err:10.4f}{abs(err / EXACT) * 100:10.1f}{-d / C[I]:14.3f}")
    return errors


def noise_simulation(sigma=0.005, trials=10000, seed=309):
    """Add Gaussian noise to every reading, recompute, take the spread."""
    random.seed(seed)
    out = {"forward": [], "backward": [], "central": []}
    for _ in range(trials):
        noisy = [c + random.gauss(0, sigma) for c in C]
        out["forward"].append(forward(noisy, I))
        out["backward"].append(backward(noisy, I))
        out["central"].append(central(noisy, I))
    print(f"\nNoise simulation: sigma = {sigma} mol/L, {trials} trials, h = 1")
    for name, vals in out.items():
        print(f"{name:10s} mean {statistics.mean(vals):8.4f}   std {statistics.stdev(vals):.4f}")
    print(f"theory:    one-sided {sigma * math.sqrt(2):.4f}, central {sigma * math.sqrt(2) / 2:.4f}")


if __name__ == "__main__":
    print("t (min):", T)
    print("C (mol/L):", C)
    e1 = report(1)
    e2 = report(2)
    print("\nError ratio (h=2 / h=1):",
          {k: round(abs(e2[k] / e1[k]), 1) for k in e1})
    noise_simulation()
