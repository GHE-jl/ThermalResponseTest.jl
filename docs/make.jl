using ThermalResponseTest
using Documenter

# Make `using ThermalResponseTest` available to every doctest in docstrings and pages.
DocMeta.setdocmeta!(
    ThermalResponseTest,
    :DocTestSetup,
    :(using ThermalResponseTest);
    recursive = true,
)

makedocs(;
    modules = [ThermalResponseTest],
    authors = "Gabriel-Dion <dion.gabriel100@gmail.com>",
    sitename = "ThermalResponseTest.jl",
    format = Documenter.HTML(;
        canonical = "https://GHE-jl.github.io/ThermalResponseTest.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
        sidebar_sitename = false,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Data & utilities" => "data_utilities.md",
        "Interpretation theory" => [
            "Overview" => "theory/overview.md",
            "First-order approximation" => "theory/first_order_approximation.md",
            "Model inversion" => "theory/model_inversion.md",
        ],
        "API reference" => "api.md",
        "References" => "references.md",
    ],
    # Keep the build strict so broken cross-references or missing docstrings fail CI.
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/GHE-jl/ThermalResponseTest.jl",
    devbranch = "main",
)
