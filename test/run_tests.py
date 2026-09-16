#!/usr/bin/env python3
"""Bounded regression and convergence checks; no output is written into the repository."""

import math
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def run(command, cwd, *, failure=False):
    result = subprocess.run(
        [str(item) for item in command], cwd=cwd, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120,
    )
    if failure:
        if result.returncode == 0 or "SolverSolveUMFPACK2 failed" not in result.stdout:
            raise RuntimeError(f"Expected solver failure: {command}\n{result.stdout}")
    elif result.returncode != 0:
        raise RuntimeError(f"Command failed: {command}\n{result.stdout[-12000:]}")
    return result.stdout


def values(output, label):
    matches = re.findall(re.escape(label) + r"\s*=?\s*(" + NUMBER + r")", output)
    if not matches:
        raise AssertionError(f"Missing {label!r}\n{output[-4000:]}")
    result = [float(item.replace("D", "E").replace("d", "e")) for item in matches]
    if not all(math.isfinite(item) for item in result):
        raise AssertionError(f"Nonfinite {label}: {result}")
    return result


def convergence(errors, minimum, label):
    if len(errors) != 3 or not all(item > 0 for item in errors):
        raise AssertionError(f"Invalid error sequence for {label}: {errors}")
    rates = [math.log2(a / b) for a, b in zip(errors, errors[1:])]
    if min(rates) < minimum:
        raise AssertionError(f"Insufficient convergence for {label}: {errors}, rates={rates}")
    print(f"PASS {label}: rates " + ", ".join(f"{rate:.3f}" for rate in rates))


def main():
    with tempfile.TemporaryDirectory(prefix="ffem-check-") as directory:
        work = Path(directory)
        (work / "output/static").mkdir(parents=True)
        output = run([ROOT / "test/test_regression"], work)
        if "PASS:" not in output:
            raise AssertionError(output)
        print("PASS core regression assertions")
        run([ROOT / "test/test_regression", "singular-fatal"], work, failure=True)
        print("PASS failed solve terminates without an explicit status argument")
        for name in ("test_mesh", "test_matvec", "test_quicksort", "test_mpi"):
            run([ROOT / "test" / name], work)
            print(f"PASS {name} smoke")

        poisson_l2, poisson_h1 = [], []
        stokes_l2, stokes_h1, pressure_l2 = [], [], []
        for n in (4, 8, 16):
            output = run([ROOT / "test/test_poisson", n, n], work)
            poisson_l2.append(values(output, "L2 error")[0])
            poisson_h1.append(values(output, "H1 error")[0])
            output = run([ROOT / "test/test_stokes", n, n], work)
            stokes_l2.append(values(output, "L2 error of u")[0])
            stokes_h1.append(values(output, "H1 error of u")[0])
            pressure_l2.append(values(output, "L2 error of p")[0])
        convergence(poisson_l2, 2.8, "Q2 Poisson L2")
        convergence(poisson_h1, 1.8, "Q2 Poisson H1")
        convergence(stokes_l2, 2.8, "P2/P1 Stokes velocity L2")
        convergence(stokes_h1, 1.8, "P2/P1 Stokes velocity H1")
        convergence(pressure_l2, 1.8, "P2/P1 Stokes pressure L2")
        for typ, order in ((121, 0.9), (122, 1.8)):
            output = run([ROOT / "test/test_interpolation", 3, typ], work)
            convergence(values(output, "Error L2"), order, f"RT {typ} L2")
            convergence(values(output, "Error Hdiv"), order, f"RT {typ} H(div)")

        for name in ("hw1", "hw2", "final"):
            run(["make", "-C", ROOT / "apps/homework" / name, "all"], work)
        run([ROOT / "apps/homework/hw1/test_jacobi_GS", 4, 4, 1, 20, "output/", 1], work)
        for name in ("hw2", "final"):
            for smoother in (1, 2):
                args = [2, 2, 4, 2, 2, smoother, f"{name}-{smoother}.txt"]
                if name == "final":
                    args += [1.0, 0.1, 1.0]
                output = run([ROOT / "apps/homework" / name / "test_multigrid", *args], work)
                residuals = values(output, "Residual =")
                if residuals[-1] >= 1e-5 * residuals[0]:
                    raise AssertionError(f"Multigrid residual did not decrease: {name}, {smoother}")
                print(f"PASS {name} multigrid smoother {smoother}")

        run(["make", "-C", ROOT / "apps/magnetic", "all"], work)
        for name in ("test_moving_exact", "test_static"):
            output = run([ROOT / "apps/magnetic" / name, 4, 0.1], work)
            values(output, "L2 error of B")
            if re.search(r"\b(?:nan|infinity)\b", output, re.IGNORECASE):
                raise AssertionError(f"Nonfinite magnetic output: {name}")
            print(f"PASS {name} smoke (not a conservation/convergence certification)")
    print("All bounded regression checks passed.")


if __name__ == "__main__":
    main()
