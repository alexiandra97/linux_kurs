# Linux config

## 1. Install tools

```bash
# install preprequisites
sudo dnf config-manager --set-enabled crb sudo dnf install git
epel-release npm
sudo dnf install ripgrep fd-find

# install nvim and tmux
sudo dnf install neovim tmux

# install oh-my-bash
bash -c "$(curl -fsSL https://raw.githubusercontent.com/ohmybash/oh-my-bash/master/tools/install.sh)"
```


