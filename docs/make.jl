using Pkg
Pkg.activate(@__DIR__)

function remove_jupyter_files()
  dirs = [joinpath(r, d) for (r, ds, _) in walkdir(@__DIR__) for d in ds if d == ".ipynb_checkpoints"]
  foreach(p -> rm(p; recursive = true), dirs)
end

remove_jupyter_files()

using Documenter
using GBGKnudsenGas

makedocs(
  sitename = "GBGKnudsenGas.jl",
  modules  = [GBGKnudsenGas],
  format   = Documenter.HTML(),
  pages    = [
    "Home" => "index.md",
    "API"  => "api.md"
  ],
)

deploydocs(
  repo      = "github.com/GrigorySarnitsky/GBGKnudsenGas.jl",
  devbranch = "master",
)
