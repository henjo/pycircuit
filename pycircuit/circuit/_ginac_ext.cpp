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
#include <pybind11/complex.h>
#include <ginac/ginac.h>
#include <complex>
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

// Compile a (complex-valued) symbolic expression to native code and evaluate it
// over a batch of real parameter points.  ``expr`` is a GiNaC-parseable string
// that may contain the imaginary unit ``I``; ``names`` are the real parameters
// (row order of ``points``).  The expression's real and imaginary parts are
// compiled together (compile_ex, native machine code) and evaluated per point.
static std::vector<std::complex<double>>
eval_sweep(const std::string &expr,
           const std::vector<std::string> &names,
           const std::vector<std::vector<double>> &points) {
    symtab tab;
    tab["I"] = I;  // imaginary unit, not a symbol
    lst syms;
    for (const auto &nm : names) {
        realsymbol sy(nm);
        tab[nm] = sy;
        syms.append(sy);
    }
    parser reader(tab);
    ex e = reader(expr);

    lst outs = {e.real_part(), e.imag_part()};
    FUNCP_CUBA fp;
    compile_ex(outs, syms, fp);  // native compile (needs a C++ compiler at runtime)

    int ndim = static_cast<int>(names.size()), ncomp = 2;
    std::vector<std::complex<double>> res;
    res.reserve(points.size());
    double out[2];
    for (const auto &pt : points) {
        fp(&ndim, pt.data(), &ncomp, out);
        res.emplace_back(out[0], out[1]);
    }
    return res;
}

PYBIND11_MODULE(_ginac_ext, m) {
    m.doc() = "Fraction-free symbolic linear solve + fast eval via GiNaC.";
    m.def("solve_numden", &solve_numden,
          py::arg("n"), py::arg("entries"), py::arg("rhs"), py::arg("symbols"));
    m.def("eval_sweep", &eval_sweep,
          py::arg("expr"), py::arg("names"), py::arg("points"));
}
