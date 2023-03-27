using MeasurementErrorModels
using Documenter

DocMeta.setdocmeta!(MeasurementErrorModels, :DocTestSetup, :(using MeasurementErrorModels); recursive=true)

makedocs(;
    modules=[MeasurementErrorModels],
    authors="Mattriks <Mattriks@users.noreply.github.com> and contributors",
    repo="https://github.com/Mattriks/MeasurementErrorModels.jl/blob/{commit}{path}#{line}",
    sitename="MeasurementErrorModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
        mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/physics"]), :tex => Dict("inlineMath" => [["\$","\$"], ["\\(","\\)"]],
        "tags" => "ams", "packages" => ["base", "ams", "autoload", "physics"],),
        )),),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Guide" => ["man/Example1.md", "man/Example2.md", "man/Example3.md"],
        "Library" => "lib/library.md"
    ],
)


deploydocs(;
    repo="github.com/Mattriks/MeasurementErrorModels.jl",
    devbranch="main",
)