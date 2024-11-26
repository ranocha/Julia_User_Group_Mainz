# VS-Code and julia over SSH

This is a tutorial on setting up usage of `julia` on a remote server using SSH via VS Code to be able to interactively execute numerical codes and visualize results of the simulations. Here are the steps we are going to take

## Install VS Code and remote SSH extension (to be run locally)

VS Code can be installed from the following link https://code.visualstudio.com/. Once you have installed VS Code, you need to install the [Remote - SSH](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh) extension by searching for its name in the extension marketplace. Within VS code, you can use the hotkey `ctrl` + `shift` + `x` for Linux/Windows or `cmd` + `shift` + `x` for `macos`. You can also follow this [guide for VS Code extensions](https://code.visualstudio.com/docs/editor/extension-marketplace) if you prefer to use GUI.

## Installing `julia` (to be done on server, NOT needed for apollo)

We need to install `julia` on the remote server. The installation is simple by using [juliaup](https://github.com/JuliaLang/juliaup) which is now the [officially recommended way](https://julialang.org/downloads/). You may need to close and reopen the shell. The installation is successful if you can start `julia` by entering `julia` in the shell. This is not needed for apollo because it has a common `julia` installed for all users.

## Connecting to remote server using SSH on VS code

The later steps can be done using GUI, but we perform it using the VS Code command palatte as it is a useful tool to learn. After opening VS code, open the VS Code command palatte by pressing `ctrl`/`cmd` (for Linux/MacOS) + `shift` + `p` for opening the command palatte. You can enter full commands here, or choose from the menu. The fastest way is to write a part of the command. In this case, enter `SSH` and choose the command `remote-SSH: Add New Host...`. Now enter the command that you use to `ssh` into the remote server. For this tutorial, we all have access to `apollo` with username `julia24`. Thus, you can enter `ssh julia24@apollo.geo.uni-mainz.de` and press `Enter`. Press `Enter` again to select the default config file (Feel free to change it if you know what you are doing). You will be given the password to `apollo` during the workshop. Once you have the password, [you are generally recommended to use a password-less login using SSH keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) although we can manage without it for the tutorial. To connect to `apollo` via SSH on VS Code, open the command palatte again (pressing `ctrl`/`cmd` (for Linux/MacOS) + `shift` + `p`) and enter/choose the command `Remote-SSH: connect current window to host` (uses current window) or `Remote-SSH: connect to host` (opens a new window). Now enter/choose `apollo.geo.uni-mainz.de`. You will need to enter the password if you haven't set up the SSH keys. Çhoose `linux` as the OS of the server if you are asked.

## Installing the `julia` extension on VS Code (to be done on server)

Once you are connected to `apollo` via SSH on VS code (you should see `SSH:apollo.geo.uni-mainz.de` at the bottom left), you should install the [`julia` extension](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia) by searching for it the extension marketplace (note that you need to do this while remaining connected to the server in VS code). We now test that the extension works by the code we will be using for this tutorial. You can download the code by using your own terminal, but we will use the VS code integrated terminal. If you don't see the integrated terminal at the bottom, you can toggle its visibility by pressing `ctrl`/`cmd` (Linux/Mac) + `j`. Then, enter the following command in the terminal (please replace `unique_name` with a unique name as all of us will be using the same server for the tutorial)
```shell
git clone https://github.com/ranocha/Julia_User_Group_Mainz.git unique_name
```
Now we will open the directory in VS code by pressing `ctrl`/`cmd` (depending on OS) + `o` and then choosing the following directory (Don't forget to replace `unique_name` with the one you chose)
```
/local/home/julia24/unique_name/Julia_User_Group_Mainz/2024-11-28__VSCode_SSH
```
If you are using some other username, replace `julia24` with that username and its location before it. Now we use `julia` through the VS code extension by again using the command palatte (`ctrl`/`cmd` + `shift` + `p`) and entering `Julia: Start REPL`. You should always make sure that you are in the current 'environment` by entering
```julia
julia> import Pkg
julia> Pkg.activate(".")
```
We now install the required dependencies for the packages by the following command. It takes some time, but is only needed once
```julia
julia> # Pkg.instantiate() (Skipped for the tutorial. Don't run this!)
```
For the tutorial, we skip this command and instead use the following command to activate an environment where I have already installed the packaged
```julia
julia> Pkg.activate("/local/home/julia24/arpit/2024-11-28__VSCode_SSH")
```
We can now run the codes in this directory. You can open the codes "my_run.jl", "elixir_advection_basic.jl", "run_blast.jl" by using file explorer. You can access the file explorer with the hotkey `ctrl`/`cmd` + `shift` + `e` or by choosing the top icon in the left side bar. Once you open these files in VS code, you can run them by pressing the play button at top right of the window, or by using the hotkeys below

## Useful VS Code hotkeys

1. `shift` + `enter` = Run `julia` code line and move to the next
2. `ctrl` + `enter` = Run `julia` code without moving to the next line
3. `alt` + `shift` + `enter` = Run block of julia code and move cursor to the next block
4. `ctrl`/`cmd` + `p` hotkey = Navigate through open folders
5. `ctrl`/`cmd` + `` ` `` = Change focus to shell
6. `ctrl`/`cmd` + 1/2/3 = change focus to active window number 1/2/3
7. `ctrl`/`cmd` + `shift` + `f` = global search in the active folder
