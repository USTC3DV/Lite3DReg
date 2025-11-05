#!/usr/bin/env bash
set -e

# 安装系统依赖
apt-get update
xargs -a apt.txt apt-get install -y

# 清理并创建 build 目录
rm -rf build
mkdir build
cd build

# 生成 Makefile，指定 Python
cmake .. -DPYTHON_EXECUTABLE=$(which python3)

# 编译
cmake --build . --config Release -j$(nproc)

# 回到主目录
cd ..

echo "✅ Build complete. .so is in build/python/"
