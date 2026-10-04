"""Device classes as wide as the C buffers hold, and one wider (stage 5,
testing for development, 2026-10-04); generated once, the signatures
written out (a class's terminals are its `analog`'s parameters).

The evaluate core's per-element buffers hold 64 unknowns (`_tran_core`:
`X[64]`, `o[4096]`; wider is `Unservable`), the walk's 64 (`_hdl_climit`,
`len(nm) <= 64`), the limiter kernel's write-back 32 (`_WB_MAX`; wider
refuses the kernel).  `Wide64`/`Wide65` and `LimWide32`/`LimWide33` sit on
each side of those limits: a ring of `k` terminals, each neighbouring pair
a cubic conductance (and, in `LimWide*`, a `limit_vds` on its voltage).
Their source is a file so the compile cache keys them (an `exec`'d class
is recompiled by every process).
"""
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, limit_vds, var
from pycircuit.utilities.param import Parameter

G = [Parameter(name='g', desc='conductance scale', unit='S', default=1e-3)]


def _ring(ts, g, limited):
    out = []
    k = len(ts)
    for j in range(k):
        b = Branch(ts[j], ts[(j + 1) % k])
        v = limit_vds(b.V) if limited else b.V
        u = var(v, f'u{j}')
        out.append(Contribution(b.I, g * (u + u * u * u)))
    return tuple(out)


class Wide64(Behavioural):
    instparams = G

    @staticmethod
    def analog(t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63):
        return _ring((t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63,), g, False)  # noqa: F821


class Wide65(Behavioural):
    instparams = G

    @staticmethod
    def analog(t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64):
        return _ring((t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62, t63, t64,), g, False)  # noqa: F821


class LimWide32(Behavioural):
    instparams = G

    @staticmethod
    def analog(t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31):
        return _ring((t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31,), g, True)  # noqa: F821


class LimWide33(Behavioural):
    instparams = G

    @staticmethod
    def analog(t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32):
        return _ring((t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32,), g, True)  # noqa: F821
