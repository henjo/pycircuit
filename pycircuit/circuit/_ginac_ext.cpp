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

// Coefficients of ``poly`` as a polynomial in ``var``, low order first
// (index k = coefficient of var^k).  ``poly`` must be expanded in ``var``.
static std::vector<std::string> coeff_list(const ex &poly, const ex &var) {
    int d = poly.degree(var);
    std::vector<std::string> c;
    c.reserve(d + 1);
    for (int k = 0; k <= d; ++k)
        c.push_back(ex_to_str(poly.coeff(var, k)));
    return c;
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

// Native solve result: keeps the numerators and the shared denominator as GiNaC
// ``ex`` objects rather than stringifying them.  For a large symbolic solve the
// determinant-sized ``num[]``/``den`` are expensive to round-trip through sympy
// (~80% of end-to-end time); callers usually want only a *compact* piece -- a
// transfer function ``x_i/x_j = num_i/num_j`` (the shared denominator cancels)
// or the denominator polynomial for poles.  Those cancellations happen here, in
// GiNaC, and only the small result crosses back to Python as a string.
struct GinacResult {
    std::vector<ex> num;  // numerators; solution component x_i = num[i] / den
    ex den;               // shared denominator (the network determinant)
    symtab tab;           // name -> symbol, to look up variables by identity

    size_t size() const { return num.size(); }

    // Raw (uncancelled) numerator / denominator -- may be large.
    std::string num_str(size_t i) const { return ex_to_str(num.at(i)); }
    std::string den_str() const { return ex_to_str(den); }

    // Full solution component x_i = num_i / den (can be large).
    std::string component_str(size_t i) const {
        return ex_to_str((num.at(i) / den).normal());
    }

    // Transfer function x_i / x_j = num_i / num_j: the shared denominator
    // cancels, so this stays compact even when the full solution is huge.
    std::string tf_str(size_t i, size_t j) const {
        return ex_to_str((num.at(i) / num.at(j)).normal());
    }

    const ex &var_ex(const std::string &var) const {
        auto it = tab.find(var);
        if (it == tab.end())
            throw std::runtime_error("unknown variable: " + var);
        return it->second;
    }

    // Transfer function x_i/x_j as polynomial-in-``var`` coefficient vectors
    // (numerator, denominator), low order first -- the canonical N(var)/D(var)
    // form.  Cancellation + collection run natively; only the (small)
    // coefficients cross back to Python.
    std::pair<std::vector<std::string>, std::vector<std::string>>
    tf_coeffs(size_t i, size_t j, const std::string &var) const {
        ex e = (num.at(i) / num.at(j)).normal();
        const ex &sv = var_ex(var);
        return {coeff_list(e.numer().expand(), sv),
                coeff_list(e.denom().expand(), sv)};
    }

    // Coefficients of the denominator (network determinant) as a polynomial in
    // ``var`` -- the characteristic polynomial whose roots are the poles.
    std::vector<std::string> denominator_coeffs(const std::string &var) const {
        return coeff_list(den.expand(), var_ex(var));
    }
};

// Solve A x = b like solve_numden, but return a native GinacResult holding the
// GiNaC ``ex`` numerators/denominator instead of stringified output.
static GinacResult
solve_native(int n,
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

    GinacResult res;
    res.tab = tab;  // keep the symbol table for later coeff/series lookups
    res.num.reserve(n);
    for (int i = 0; i < n; ++i)
        res.num.push_back((x(i, 0) * det).normal());
    res.den = det.normal();
    return res;
}

// Polynomial-in-``var`` coefficients of an arbitrary rational expression
// ``expr``: returns (numerator_coeffs, denominator_coeffs), low order first.
// normal()/numer_denom()/collection run in GiNaC; only coefficients are
// stringified.  ``symbols`` lists every symbol in ``expr`` (incl. ``var``).
static std::pair<std::vector<std::string>, std::vector<std::string>>
poly_coeffs(const std::string &expr,
            const std::string &var,
            const std::vector<std::string> &symbols) {
    symtab tab;
    for (const auto &nm : symbols)
        tab[nm] = symbol(nm);
    if (!tab.count(var))
        tab[var] = symbol(var);
    parser reader(tab);
    ex e = reader(expr).normal();
    ex sv = tab[var];
    return {coeff_list(e.numer().expand(), sv),
            coeff_list(e.denom().expand(), sv)};
}

// Truncated series of ``expr`` in ``var`` about ``point`` to ``order`` terms,
// returned as a single ordinary-polynomial string (the Order term dropped).
// Laurent expansion for low-/high-frequency and pole/zero approximation.
static std::string
series_expr(const std::string &expr,
            const std::string &var,
            const std::string &point,
            int order,
            const std::vector<std::string> &symbols) {
    symtab tab;
    for (const auto &nm : symbols)
        tab[nm] = symbol(nm);
    if (!tab.count(var))
        tab[var] = symbol(var);
    parser reader(tab);
    ex e = reader(expr);
    ex sv = tab[var];
    ex pt = reader(point);
    ex se = series_to_poly(e.series(sv == pt, order));
    return ex_to_str(se);
}

PYBIND11_MODULE(_ginac_ext, m) {
    m.doc() = "Fraction-free symbolic linear solve + fast eval via GiNaC.";
    m.def("solve_numden", &solve_numden,
          py::arg("n"), py::arg("entries"), py::arg("rhs"), py::arg("symbols"));
    m.def("eval_sweep", &eval_sweep,
          py::arg("expr"), py::arg("names"), py::arg("points"));

    py::class_<GinacResult>(m, "GinacResult",
          "Native GiNaC solve result (num[]/den kept as GiNaC ex).")
        .def("__len__", &GinacResult::size)
        .def("num_str", &GinacResult::num_str, py::arg("i"))
        .def("den_str", &GinacResult::den_str)
        .def("component_str", &GinacResult::component_str, py::arg("i"))
        .def("tf_str", &GinacResult::tf_str, py::arg("i"), py::arg("j"))
        .def("tf_coeffs", &GinacResult::tf_coeffs,
             py::arg("i"), py::arg("j"), py::arg("var"))
        .def("denominator_coeffs", &GinacResult::denominator_coeffs,
             py::arg("var"));
    m.def("solve_native", &solve_native,
          py::arg("n"), py::arg("entries"), py::arg("rhs"), py::arg("symbols"));
    m.def("poly_coeffs", &poly_coeffs,
          py::arg("expr"), py::arg("var"), py::arg("symbols"));
    m.def("series_expr", &series_expr,
          py::arg("expr"), py::arg("var"), py::arg("point"), py::arg("order"),
          py::arg("symbols"));
}
