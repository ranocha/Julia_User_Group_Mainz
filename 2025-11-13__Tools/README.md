# Tools for programming and related tasks

This continues the discussions we started in the summer term 2025 on
Julia workflows [here](https://github.com/ranocha/Julia_User_Group_Mainz/tree/main/2025-05-08__Julia_workflows).

For teaching, I still use mainly [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebooks. Nowadays, I have decided to share them via GitHub instead of Moodle, e.g.,
- [Mathematics for Computer Scientists 1, Winter Term 2025/2026](https://github.com/ranocha/2025_MfI1)
- [Mathematics for Computer Scientists 2b, Winter Term 2024/2025](https://github.com/ranocha/2024_MfI2b)
- [Introduction to Numerical Mathematics, Summer Term 2025](https://github.com/ranocha/2025_Num1)
- [Numerical Methods for Ordinary Differential Equations, Winter Term 2025/2026](https://github.com/ranocha/2025_NumODE)

When developing code (for papers), I write Julia scripts using a setup like

I typically use [Visual Studio Code](https://code.visualstudio.com/)
to develop the code and write the paper (using the LaTeX and Julia
extensions). In addition, I often use a standard terminal next to
VSCode where I start a REPL session as

```julia
julia> using Revise; includet("code/code.jl")
```

Then, I can modify the source code and use the latest version in the
Julia REPL. The `code.jl` file typically begins with something like

```julia
# Install packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Load packages
using Trixi
```

When working on packages, I prefer a similar workflow - VSCode plus
terminal. Inside the package directory, I create a `run` directory.
While working on the package, I start Julia with the `run` directory
as project where I `Pkg.develop` the package I am working on. This
allows me to `Pkg.add` additional packages I use for development.


## Tools

Here is a list of tools I find useful. Please feel free to suggest
others via PRs to this repository.


### Julia Extension for VSCode

The [Julia extension for VSCode](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia)
is very useful for writing and debugging Julia code.


### LaTeX Workshop

I use the [LaTeX Workshop](https://marketplace.visualstudio.com/items?itemName=James-Yu.latex-workshop)
extension with VSCode for writing LaTeX documents
such as papers or lecture notes.


### LTeX – LanguageTool grammar/spell checking

The [LTeX extension](https://marketplace.visualstudio.com/items?itemName=valentjn.vscode-ltex)
for VSCode provides grammar and spell checking for
various documents including LaTeX documents, Markdown
files etc.

Let's find typos in the notes from 2025-05-08.


### Spell Checker CI: typos

The [typos](https://github.com/crate-ci/typos)
spell checker provides a GitHub action that we use
in various repositories, e.g.,
<https://github.com/trixi-framework/Trixi.jl/blob/4e3926430471902ba7e541f258e6f5324ea1fa54/.github/workflows/SpellCheck.yml>.
You can also install it locally, e.g., via [homebrew](https://formulae.brew.sh/formula/typos-cli).


### GitHub Copilot

We have seen [GitHub copilot](https://github.com/features/copilot)
in action, e.g., in the presentation of Evangelos Moulas some time ago.
It is accessible via a VSCode extension.

GitHub made it accessible for students/teachers via its education program
[already some time ago](https://github.blog/news-insights/product-news/github-copilot-now-available-for-teachers/).
It can not only complete code, but also LaTeX documents, e.g.,
lecture notes.


### Coding Agents

GitHub also provides access to various coding agents, e.g., Claude.
They can be useful for tedious tasks such as writing extensive tests.
I used Claude successfully while working on the issue
<https://github.com/ranocha/BSeries.jl/issues/285>.
The related PRs are
- <https://github.com/ranocha/BSeries.jl/pull/286>
- <https://github.com/ranocha/BSeries.jl/pull/287>

I first wrote the code with examples in the docstrings that I checked
carefully. Then, I asked Claude to generate tests and checked them.
This saved me quite some time.

To open the chat with a coding agent, you can hit `Shift + Cmd + I` (macOS) in VSCode.

Chris Rackauckas described his use of Claude in

- <https://discourse.julialang.org/t/the-use-of-claude-code-in-sciml-repos/131009/8>
- <https://www.stochasticlifestyle.com/claude-code-in-scientific-computing-experiences-maintaining-julias-sciml-infrastructure/>


### ChatGPT

[ChatGPT](https://chatgpt.com) has been around for some time and I assume you all know about it.


### Your experiences?
