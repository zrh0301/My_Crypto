# /home/zrheng/Documents/My_Crypto/sm2/Makefile

# ==========================================
# 编译器设置
# ==========================================
CC = gcc
CXX = g++

# ==========================================
# 路径设置 (Directory Definitions)
# ==========================================
# BN (Big Number) 模块路径
BN_INC_DIR = ../bn/inc
BN_SRC_DIR = ../bn/src

# ECC (Elliptic Curve) 模块路径
# [新增] 根据你的描述添加了 ECC 的路径
ECC_INC_DIR = ../ecc/inc
ECC_SRC_DIR = ../ecc/src

# SM2 (当前目录) 路径
# 当前目录直接用 . 即可

# ==========================================
# 编译选项 (Flags)
# ==========================================
# [关键] -I 选项告诉编译器去哪里找 .h 头文件
# 我们需要包含 BN 和 ECC 的头文件目录
INCLUDES = -I$(BN_INC_DIR) -I$(ECC_INC_DIR) -I.

CFLAGS   = $(INCLUDES) -Wall -O3
CXXFLAGS = $(INCLUDES) -Wall -O3

# ==========================================
# 目标与依赖 (Targets & Objects)
# ==========================================

# 最终生成的可执行文件
TARGET = sm2_test

# 对象文件列表
# 我们需要：bn 的 C 对象，ecc 的 C++ 对象，以及 sm2 自己的对象
OBJS = wlc_bn.o ecc.o sm2.o

# ==========================================
# 构建规则 (Build Rules)
# ==========================================

# 默认目标
all: $(TARGET)

# 1. 链接 (Linking)
# 将所有生成的 .o 文件链接成可执行文件
$(TARGET): $(OBJS)
	$(CXX) $^ -o $@

# 2. 编译 wlc_bn.c (位于 ../bn/src)
# 依赖: wlc_bn.c 和 wlc_bn.h
wlc_bn.o: $(BN_SRC_DIR)/wlc_bn.c $(BN_INC_DIR)/wlc_bn.h
	$(CC) $(CFLAGS) -c $< -o $@

# 3. 编译 ecc.cpp (位于 ../ecc/src)
# [新增规则] 告诉 make 去哪里找 ecc.cpp
# 依赖: ecc.cpp, ecc.h, 以及底层的 wlc_bn.h
ecc.o: $(ECC_SRC_DIR)/ecc.cpp $(ECC_INC_DIR)/ecc.h $(BN_INC_DIR)/wlc_bn.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 4. 编译 sm2.cpp (位于当前目录)
# 依赖: sm2.cpp, ecc.h, wlc_bn.h
sm2.o: sm2.cpp $(ECC_INC_DIR)/ecc.h $(BN_INC_DIR)/wlc_bn.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==========================================
# 清理与重构
# ==========================================
clean:
	rm -f $(OBJS) $(TARGET)

rebuild: clean all

.PHONY: all clean rebuild