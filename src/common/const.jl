# using Unitful
# using UnitfulAtomic

# const Bohr = uconvert(u"Å", 1u"bohr").val
# const kB = uconvert(u"eV/K", 1u"k").val

using WannierIO: Bohr
const kB = 8.617333262145179e-5  # eV/K
