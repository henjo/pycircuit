// pybind11 extension: fraction-free symbolic linear solve via GiNaC.
//
// Exposes one function, solve_numden(n, entries, rhs, symbols), used by
// pycircuit.circuit._ginac / GinacToolkit.  All GiNaC contact stays in C++;
// the Python side passes/receives GiNaC-parseable strings.  The win over
// sympy is GiNaC's normal(): it returns compact, GCD-cancelled polynomials.
//
// Build (needs GiNaC via pkg-config + pybind11 + python headers):
//   g++ -O2 -std=c++17 -shared -fPIC $(python3-config --includes) \
//       -I<pybind11 include> _ginac_ext.cpp -o _ginac_ext.so \
//       $(pkg-config --cflags --libs ginac)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <ginac/ginac.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;
using namespace GiNaC;

static std::string ex_to_str(const ex &e) {
    std::ostringstream os;
    os << e;
    return os.str();
}

// Solve A x = b for an n x n system whose entries are GiNaC-parseable strings.
// Returns (numerators, denominator) as strings with the SHARED denominator equal
// to the network determinant, so x_i = numerators[i] / denominator.  Raises on a
// parse error or a singular system; the Python side falls back to sympy.
static std::pair<std::vector<std::string>, std::string>
solve_numden(int n,
             const std::vector<std::string> &entries,  // n*n, row-major
             const std::vector<std::string> &rhs,       // n
             const std::vector<std::string> &symbols) {
    symtab tab;
    for (const auto &nm : symbols)
        tab[nm] = symbol(nm);
    parser reader(tab);

    matrix A(n, n), b(n, 1), vars(n, 1);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A.set(i, j, reader(entries[i * n + j]));
    for (int i = 0; i < n; ++i) {
        b.set(i, 0, reader(rhs[i]));
        vars.set(i, 0, symbol());
    }

    ex det = A.determinant();
    matrix x = A.solve(vars, b);  // throws std::runtime_error if singular

    std::vector<std::string> nums;
    nums.reserve(n);
    for (int i = 0; i < n; ++i)
        nums.push_back(ex_to_str((x(i, 0) * det).normal()));
    return {nums, ex_to_str(det.normal())};
}

PYBIND11_MODULE(_ginac_ext, m) {
    m.doc() = "Fraction-free symbolic linear solve via GiNaC (normal()).";
    m.def("solve_numden", &solve_numden,
          py::arg("n"), py::arg("entries"), py::arg("rhs"), py::arg("symbols"));
}
