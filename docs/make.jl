using Pkg
Pkg.activate(@__DIR__)

function remove_jupyter_files()
  dirs = [joinpath(r, d) for (r, ds, _) in walkdir(@__DIR__) for d in ds if d == ".ipynb_checkpoints"]
  foreach(p -> rm(p; recursive = true), dirs)
end

remove_jupyter_files()

using Documenter
using DocumenterMarkdown
using GBGKnudsenGas

makedocs(
  sitename = "GBGKnudsenGas.jl",
  modules  = [GBGKnudsenGas],
  pages    = ["Home" => "index.md"],
  format   = Markdown(),
  #warnonly = true,   # optional: turn errors into warnings while experimenting
)

cp(joinpath(@__DIR__, "build", "index.md"), joinpath(@__DIR__, "..", "README.md"); force = true)