# Build the relativistic (Teukolsky) GW correction tables, one channel per invocation:
#   julia scripts/build_gw_kerr_tables.jl --kind tr --a 3-2-2 --b 2-1-1
#   julia scripts/build_gw_kerr_tables.jl --default --lmax_levels 5 --task 7 --ntasks 60
#   julia scripts/build_gw_kerr_tables.jl --default --list          (print the default channel set)
# Writes rate_sve/gw_kerr/<kind>_<a>_<b>.dat (read by gw_build_cache with gw_model=:rel).
# Grid: the superradiance (α, a) grid subsampled (--alpha_stride, --spin_stride), α >= --alpha_min.
using ArgParse, Suppressor
include(joinpath(@__DIR__, "..", "Core", "constants.jl"))
include(joinpath(@__DIR__, "..", "state_utils.jl"))
@suppress include(joinpath(@__DIR__, "..", "Numerics", "gw_rates.jl"))
include(joinpath(@__DIR__, "..", "Numerics", "gw_kerr_build.jl"))

s = ArgParseSettings()
@add_arg_table! s begin
    "--kind";          arg_type = String;  default = "";   help = "tr | ann"
    "--a";             arg_type = String;  default = "";   help = "level a as n-l-m"
    "--b";             arg_type = String;  default = "";   help = "level b as n-l-m"
    "--default";       action = :store_true;              help = "use gw_kerr_default_channels (--lmax_levels, --overtones)"
    "--lmax_levels";   arg_type = Int;     default = 5
    "--overtones";     arg_type = Int;     default = 1
    "--task";          arg_type = Int;     default = 0;    help = "with --default: 1-based index into the channel list (0 = all)"
    "--ntasks";        arg_type = Int;     default = 1;    help = "with --default and --task: channels task, task+ntasks, ..."
    "--list";          action = :store_true;              help = "print the channel list and exit"
    "--alpha_stride";  arg_type = Int;     default = 2
    "--spin_stride";   arg_type = Int;     default = 4
    "--extra_spins";   arg_type = String;  default = "0.95"; help = "extra spin nodes, comma separated (\"\" = none); resumed tables gain the new columns"
    "--alpha_min";     arg_type = Float64; default = 0.03; help = "lowest α computed (annihilations: 0.05 recommended in Float64)"
    "--alpha_max";     arg_type = Float64; default = 0.0;  help = "override the pair's α_max (0 = superradiance-table range)"
    "--bigfloat_below"; arg_type = Float64; default = 0.0; help = "BigFloat Teukolsky integration for α below this"
    "--lmax";          arg_type = Int;     default = 0;    help = "largest ℓ̃ (0 = |m̃| + 10)"
    "--ltol";          arg_type = Float64; default = 1e-8
    "--decay";         arg_type = Float64; default = 60.0
    "--no_horizon";    action = :store_true;              help = "skip the horizon flux (Gamma_H column = NaN)"
    "--outdir";        arg_type = String;  default = "";   help = "table directory (default rate_sve/gw_kerr)"
    "--no_resume";     action = :store_true
    "--quiet";         action = :store_true
end
args = parse_args(s)

channels = Tuple{Symbol,NTuple{3,Int},NTuple{3,Int}}[]
if args["default"]
    channels = gw_kerr_default_channels(lmax=args["lmax_levels"], overtones=args["overtones"])
    if args["task"] > 0
        channels = channels[args["task"]:args["ntasks"]:end]
    end
else
    isempty(args["kind"]) && error("give --kind/--a/--b or --default")
    push!(channels, (Symbol(args["kind"]), gw_kerr_parse_level(args["a"]), gw_kerr_parse_level(args["b"])))
end
if args["list"]
    for (i, (k, sa, sb)) in enumerate(channels)
        println(i, " ", k, " ", gw_kerr_level_str(sa), " ", gw_kerr_level_str(sb))
    end
    exit(0)
end
outdir = isempty(args["outdir"]) ? gw_kerr_table_dir() : args["outdir"]
for (k, sa, sb) in channels
    gw_kerr_build_channel(k, sa, sb; path=joinpath(outdir, gw_kerr_table_name(k, sa, sb)), resume=!args["no_resume"],
                          alpha_stride=args["alpha_stride"], spin_stride=args["spin_stride"], alpha_min=args["alpha_min"],
                          extra_spins=[parse(Float64, x) for x in split(args["extra_spins"], ',') if !isempty(strip(x))],
                          alpha_max=args["alpha_max"] > 0 ? args["alpha_max"] : nothing,
                          bigfloat_below=args["bigfloat_below"], lmax=args["lmax"] > 0 ? args["lmax"] : nothing,
                          ltol=args["ltol"], decay=args["decay"], horizon=!args["no_horizon"], verbose=!args["quiet"])
end
