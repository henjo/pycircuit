"""Machine-made changes of one C line, as pure functions of its text (tests/test_gate_mutate.py): what a mutation sample draws from.

Each is the kind of slip a person makes: a comparison off by its equal, a condition negated, `&&` for `||`, a sign or an operation swapped, a constant off
by one, a return value flipped, a statement dropped.  Whether one compiles is the compiler's to say (`mutants.py` builds it first); whether one changes
anything is what the sample measures: some are equivalent, and only reading a survivor tells.
"""
import random
import re

_COMMENT = re.compile(r'/\*.*?\*/|//.*$')
_OPEN_COMMENT = re.compile(r'/\*.*$')
_CMP = {'<=': '<', '>=': '>', '==': '!=', '!=': '==', '<': '<=', '>': '>='}
_CMP_AT = re.compile(r'(?<![<>\-=!])(<=|>=|==|!=|<(?![<=])|>(?![>=]))')
_LOGIC = re.compile(r'&&|\|\|')
_ARITH = re.compile(r'(?<=[\w\)\]] )([-+*/])(?= [\w\(\-])')
_NUMBER = re.compile(r'(?<![\w.])(\d+\.\d*(?:[eE][-+]?\d+)?|\d+\.?[eE][-+]?\d+|\d+)(?![\w.])')
_RETURN = re.compile(r'\breturn ([01]);')
_DECL = re.compile(r'^\s*(const\b|static\b|struct\b|enum\b|typedef\b|unsigned\b|signed\b|double\b|float\b|long\b|int\b|char\b|void\b|size_t\b|\w+_t\b|return\b|break\b|continue\b|goto\b|case\b|default\b|else\b|if\b|for\b|while\b|do\b|switch\b)')
_STATEMENT = re.compile(r'^\s*[^;{}]*;\s*$')


def _code(line):
    """The line with its comments blanked, columns kept; None for a line a change cannot be read on (a preprocessor line, a comment's body)."""
    s = line.lstrip()
    if not s or s.startswith(('#', '*', '/*', '//')) or '"' in line:
        return None
    code = _COMMENT.sub(lambda m: ' ' * len(m.group()), line)
    return _OPEN_COMMENT.sub(lambda m: ' ' * len(m.group()), code)


def _balanced(text):
    depth = 0
    for ch in text:
        depth += (ch == '(') - (ch == ')')
        if depth < 0:
            return False
    return depth == 0


def _swap(line, code, pattern, change):
    return [line[:m.start()] + change(m) + line[m.end():] for m in pattern.finditer(code)]


def _number(m):
    t = m.group()
    if t.isdigit():
        return str(int(t) + 1)
    v = float(t)
    return repr(2.0 * v) if v else '1.0'


def _negated(line, code):
    """`if (X) ...` / `while (X) ...` with X's parenthesis closed on this line: `if (!(X)) ...`."""
    m = re.search(r'\b(if|while) \(', code)
    if not m:
        return []
    depth, i = 1, m.end()
    while i < len(code) and depth:
        depth += (code[i] == '(') - (code[i] == ')')
        i += 1
    if depth:
        return []
    return [line[:m.end()] + '!(' + line[m.end():i - 1] + '))' + line[i:]]


def mutations(line):
    """Every single change this line allows: `[(kind, new line)]`; empty where none reads (a preprocessor line, a comment, a line with a string)."""
    code = _code(line)
    if code is None:
        return []
    out = [('cmp', x) for x in _swap(line, code, _CMP_AT, lambda m: _CMP[m.group()])]
    out += [('neg', x) for x in _negated(line, code)]
    out += [('logic', x) for x in _swap(line, code, _LOGIC, lambda m: '||' if m.group() == '&&' else '&&')]
    out += [('arith', x) for x in _swap(line, code, _ARITH, lambda m: {'+': '-', '-': '+', '*': '/', '/': '*'}[m.group(1)])]
    out += [('ret', x) for x in _swap(line, code, _RETURN, lambda m: f'return {1 - int(m.group(1))};')]
    out += [('const', x) for x in _swap(line, code, _NUMBER, _number)]
    stripped = code.rstrip()
    if _STATEMENT.match(stripped) and not _DECL.match(stripped) and _balanced(stripped):
        indent = line[:len(line) - len(line.lstrip())]
        out.append(('drop', indent + ';'))
    seen = set()
    return [(k, x) for k, x in out if not (x in seen or seen.add(x))]


def draw(lines, n, seed, least=3):
    """`n` mutants, spread over the units by their executable lines (`least` at the fewest), reproducible from `seed`: `lines` is `{unit: {line number: text}}`;
    `[{'unit', 'line', 'kind', 'old', 'new'}]`, one a source line, in draw order: more than `n` are drawn per unit so that one that does not compile can
    be replaced (`reserve`), the first `n` of a unit's being the sample."""
    rng = random.Random(seed)
    total = sum(len(v) for v in lines.values())
    out = {}
    for unit in sorted(lines):
        want = max(least, round(n * len(lines[unit]) / total))
        order = sorted(lines[unit])
        rng.shuffle(order)
        got = []
        for ln in order:
            ms = mutations(lines[unit][ln])
            if ms:
                kinds = sorted({k for k, _x in ms})
                kind = rng.choice(kinds)
                new = rng.choice([x for k, x in ms if k == kind])
                got.append({'unit': unit, 'line': ln, 'kind': kind, 'old': lines[unit][ln], 'new': new})
        out[unit] = {'want': want, 'drawn': got}
    return out


#: a line's guard: an `if` whose whole statement is a return of a hand-back
_GUARD_RETURN = re.compile(r'^\s*return\s+([^;]*(?:BAIL|_bail\()[^;]*);')


def forced_guard(line, site):
    """A one-line hand-back check, `... if (X) return BAIL...;`, with `pcx_fault(site)` asked first: `if (pcx_fault(site) || (X)) return BAIL...;`.  None for any other line.
    A hand-back's contract is that the run goes on as if the C had never tried, so a build of these can fire each check where it never fires by itself."""
    code = _code(line)
    if code is None:
        return None
    m = re.search(r'\bif \(', code)
    if not m:
        return None
    depth, i = 1, m.end()
    while i < len(code) and depth:
        depth += (code[i] == '(') - (code[i] == ')')
        i += 1
    if depth or not _GUARD_RETURN.match(code[i:]):
        return None
    return f'{line[:m.end()]}pcx_fault({site}) || ({line[m.end():i - 1]})){line[i:]}'


def force_guards(text, env):
    """A unit's text with every one-line hand-back check forced by `pcx_fault`, which reads `env` (`line:nth`: the check on that line of the original, at its nth asking):
    `(text, {line: the check as it stood})`, the original's line numbers the sites.  The text unchanged, and no sites, where it holds no such check."""
    out, sites = [], {}
    for k, line in enumerate(text.split('\n'), 1):
        new = forced_guard(line, k)
        if new is not None:
            sites[k] = line
        out.append(line if new is None else new)
    if not sites:
        return text, {}
    head = ('#include <stdio.h>\n#include <stdlib.h>\n'
            'static long pcx_fault_line = -2, pcx_fault_at = 0, pcx_fault_seen = 0;\n'
            'static int pcx_fault(long line)\n{\n'
            '    if (pcx_fault_line == -2) {\n'
            f'        const char *e = getenv("{env}");\n'
            '        pcx_fault_line = -1;\n'
            '        if (e) sscanf(e, "%ld:%ld", &pcx_fault_line, &pcx_fault_at);\n'
            '    }\n'
            '    if (line != pcx_fault_line) return 0;\n'
            '    if (__atomic_add_fetch(&pcx_fault_seen, 1, __ATOMIC_RELAXED) != pcx_fault_at) return 0;\n'
            '    fputs("PCX_FAULT_FIRED\\n", stderr);\n'
            '    return 1;\n}\n')
    return head + '\n'.join(out), sites
