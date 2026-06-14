# theme.jl  — include after `using Plots`
using Colors

const ACCENT = colorant"#9B001F"   # JGUMITRED: accent, hero series, ROCK2 fill
const REF    = colorant"#565B66"   # reference / literature data

const SERIES = [ACCENT, colorant"#1F6F8B", colorant"#4C9F70",
                colorant"#C77F0A", colorant"#6A4C93", REF]

# const HEAT = cgrad([colorant"#F2C14E", colorant"#E07A3F",
#                     ACCENT, colorant"#4A0E1A"])

Plots.default(
    palette    = SERIES,
    fontfamily = "Computer Modern"
)