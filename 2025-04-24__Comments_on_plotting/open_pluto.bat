cd /d "%~dp0"
julia -e "import Pkg; Pkg.activate(pwd()); Pkg.instantiate(); import Pluto; Pluto.run()

