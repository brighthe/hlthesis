#!/bin/bash

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Neovim Configuration Installer ===${NC}"

# 获取脚本所在目录的绝对路径
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NVIM_CONFIG_SOURCE="$SCRIPT_DIR/nvim" # 配置文件在 nvim 文件夹中
NVIM_CONFIG_TARGET="$HOME/.config/nvim"

# 检查源配置是否存在
if [ ! -d "$NVIM_CONFIG_SOURCE" ]; then
  echo -e "${RED}Error: nvim directory not found in $SCRIPT_DIR${NC}"
  echo -e "${YELLOW}Expected structure: $SCRIPT_DIR/nvim/${NC}"
  exit 1
fi

# 显示将要安装的配置文件
echo -e "${YELLOW}Configuration files to install:${NC}"
ls -la "$NVIM_CONFIG_SOURCE"

# 确保 .config 目录存在
mkdir -p ~/.config

# 处理现有配置
if [ -d "$NVIM_CONFIG_TARGET" ] || [ -L "$NVIM_CONFIG_TARGET" ]; then
  BACKUP_NAME="nvim.backup.$(date +%Y%m%d_%H%M%S)"
  echo -e "${YELLOW}Existing nvim config found, backing up to ~/.config/$BACKUP_NAME${NC}"
  mv "$NVIM_CONFIG_TARGET" "$HOME/.config/$BACKUP_NAME"
fi

# 创建软链接
ln -sf "$NVIM_CONFIG_SOURCE" "$NVIM_CONFIG_TARGET"

# 验证安装
if [ -L "$NVIM_CONFIG_TARGET" ] && [ -d "$NVIM_CONFIG_TARGET" ]; then
  echo -e "${GREEN}✓ Neovim configuration installed successfully!${NC}"
  echo -e "Config location: ${YELLOW}$NVIM_CONFIG_TARGET${NC} -> ${YELLOW}$NVIM_CONFIG_SOURCE${NC}"

  # 显示安装的配置文件
  echo -e "\n${GREEN}Installed configuration files:${NC}"
  ls -la "$NVIM_CONFIG_TARGET"
else
  echo -e "${RED}✗ Installation failed${NC}"
  exit 1
fi

echo -e "\n${GREEN}Installation complete! You can now use nvim with your custom configuration.${NC}"
echo -e "${YELLOW}Note: Make sure you have all required dependencies installed (Node.js, Python, etc.)${NC}"
